using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    internal static class NurbsCurve2DHelper
    {
        public static Vector2 GetPointOnCurve(NurbsCurve<Vector2> curve, double paramT)
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            if (MathUtils.IsAlmostEqualTo(paramT, knots[0]))
            {
                return controlPoints[0].Value;
            }

            int n = controlPoints.Count - 1;
            if (MathUtils.IsAlmostEqualTo(paramT, knots[n + degree + 1]))
            {
                return controlPoints[n].Value;
            }

            int knotSpanIndex = Polynomials.GetKnotSpanIndex(degree, knots, paramT);
            int originMultiplicity = Polynomials.GetKnotMultiplicity(knots, paramT);

            int times = degree - originMultiplicity;

            // working in Homogeneous coordinates (w*x, w*y, w)
            var temp = new Vector3[times + 1];

            for (int i = 0; i <= times; i++)
            {
                var cp = controlPoints[knotSpanIndex - degree + i];
                double w = cp.Weight;
                temp[i] = new Vector3(cp.Value * (float)w, (float)w);
            }

            for (int j = 1; j <= times; j++)
            {
                for (int i = 0; i <= times - j; i++)
                {
                    double num = paramT - knots[knotSpanIndex - degree + j + i];
                    double den =
                        knots[i + knotSpanIndex + 1] - knots[knotSpanIndex - degree + j + i];
                    double alpha = num / den;

                    temp[i] = (float)alpha * temp[i + 1] + (float)(1.0 - alpha) * temp[i];
                }
            }

            Vector3 result = temp[0];
            if (MathUtils.IsZero(result.Z))
                return Vector2.Zero;

            return new Vector2(result.X, result.Y) / result.Z;
        }

        public static Vector2 GetPointOnCurveByCornerCut(NurbsCurve<Vector2> curve, double paramT)
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            if (MathUtils.IsAlmostEqualTo(paramT, knots[0]))
                return controlPoints[0].Value;

            int n = controlPoints.Count - 1;
            if (MathUtils.IsAlmostEqualTo(paramT, knots[n + degree + 1]))
                return controlPoints[n].Value;

            int knotSpanIndex = Polynomials.GetKnotSpanIndex(degree, knots, paramT);
            int originMultiplicity = Polynomials.GetKnotMultiplicity(knots, paramT);

            int times = degree - originMultiplicity;
            var temp = new Vector3[times + 1];

            for (int i = 0; i <= times; i++)
            {
                var cp = controlPoints[knotSpanIndex - degree + i];
                // Convert to Homogeneous
                temp[i] = new Vector3(cp.Value * (float)cp.Weight, (float)cp.Weight);
            }

            for (int j = 1; j <= times; j++)
            {
                for (int i = 0; i <= times - j; i++)
                {
                    double num = paramT - knots[knotSpanIndex - degree + j + i];
                    double den =
                        knots[i + knotSpanIndex + 1] - knots[knotSpanIndex - degree + j + i];
                    double alpha = num / den;

                    temp[i] = (float)alpha * temp[i + 1] + (float)(1.0 - alpha) * temp[i];
                }
            }

            var result = temp[0];
            if (MathUtils.IsZero(result.Z))
                return Vector2.Zero;
            return new Vector2(result.X, result.Y) / result.Z;
        }

        public static IReadOnlyList<Vector2> ComputeRationalCurveDerivatives(
            NurbsCurve<Vector2> curve,
            int derivative,
            double paramT
        )
        {
            // Homogeneous Coordinates for 2D are Vector3 (x*w, y*w, w)
            // Allocate a span-compatible array for the intermediate homogeneous derivatives
            var homogeneousDerivatives = new Vector3[derivative + 1];

            ComputeDerivatives(curve, derivative, paramT, homogeneousDerivatives);

            var derivatives = new List<Vector2>(derivative + 1);
            for (int i = 0; i <= derivative; i++)
                derivatives.Add(Vector2.Zero);

            for (int k = 0; k <= derivative; k++)
            {
                var v = new Vector2(homogeneousDerivatives[k].X, homogeneousDerivatives[k].Y);

                for (int i = 1; i <= k; i++)
                {
                    double binom = MathUtils.Binomial(k, i);
                    // The weight is stored in Z for 2D curve homogeneous coordinates
                    double weightDerivative = homogeneousDerivatives[i].Z;

                    v -= (float)(binom * weightDerivative) * derivatives[k - i];
                }

                // Divide by the weight (0-th derivative of weight is the weight itself)
                derivatives[k] = v / homogeneousDerivatives[0].Z;
            }

            return derivatives;
        }

        public static IReadOnlyList<BezierCurve<Vector2>> DecomposeToBeziers(
            NurbsCurve<Vector2> curve
        )
        {
            int degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            for (int i = 0; i < knots.Count; i++)
            {
                int multi = Polynomials.GetKnotMultiplicity(knots, knots[i]);
                Validate.Argument(
                    multi <= degree + 1,
                    "curve",
                    "curve knot multiplicity must be less than or equals to degree + 1."
                );
            }

            int n = controlPoints.Count - 1;
            int m = n + degree + 1;

            int breaks = 0;
            for (int i = 0; i <= m; ++i)
            {
                if (i == 0 || !MathUtils.IsAlmostEqualTo(knots[i], knots[i - 1]))
                {
                    breaks++;
                }
            }

            int bezierSize = breaks - 1;
            var beziersCP = new Vector3[bezierSize][];
            for (int i = 0; i < bezierSize; i++)
            {
                beziersCP[i] = new Vector3[degree + 1];
            }

            int a = degree;
            int b = degree + 1;
            int nb = 0;

            for (int i = 0; i <= degree; i++)
            {
                var cp = controlPoints[i];
                beziersCP[nb][i] = new Vector3(cp.Value * (float)cp.Weight, (float)cp.Weight);
            }

            while (b < m)
            {
                int i = b;
                while (b < m && MathUtils.IsAlmostEqualTo(knots[b + 1], knots[b]))
                {
                    b++;
                }
                int multi = b - i + 1;
                if (multi < degree)
                {
                    double numerator = knots[b] - knots[a];
                    var alphaVector = new double[degree + 1];
                    for (int j = degree; j > multi; j--)
                    {
                        alphaVector[j - multi - 1] = numerator / (knots[a + j] - knots[a]);
                    }

                    int r = degree - multi;
                    for (int j = 1; j <= r; j++)
                    {
                        int save = r - j;
                        int s = multi + j;
                        for (int k = degree; k >= s; k--)
                        {
                            double alpha = alphaVector[k - s];
                            beziersCP[nb][k] =
                                (float)alpha * beziersCP[nb][k]
                                + (float)(1.0 - alpha) * beziersCP[nb][k - 1];
                        }

                        if (b < m && nb + 1 < bezierSize)
                        {
                            beziersCP[nb + 1][save] = beziersCP[nb][degree];
                        }
                    }

                    nb++;
                    if (b < m && nb < bezierSize)
                    {
                        for (int k = degree - multi; k <= degree; k++)
                        {
                            var cp = controlPoints[b - degree + k];
                            beziersCP[nb][k] = new Vector3(
                                cp.Value * (float)cp.Weight,
                                (float)cp.Weight
                            );
                        }

                        a = b;
                        b += 1;
                    }
                }
            }

            var result = new List<BezierCurve<Vector2>>(bezierSize);
            for (int i = 0; i < bezierSize; i++)
            {
                var pts = new ControlPoint<Vector2>[degree + 1];
                for (int k = 0; k <= degree; k++)
                {
                    var h = beziersCP[i][k];
                    double w = h.Z;
                    if (MathUtils.IsZero(w))
                    {
                        pts[k] = new ControlPoint<Vector2>(Vector2.Zero, 0);
                    }
                    else
                    {
                        pts[k] = new ControlPoint<Vector2>(new Vector2(h.X, h.Y) / (float)w, w);
                    }
                }
                result.Add(new BezierCurve<Vector2>(pts));
            }

            return result;
        }

        public static (
            IReadOnlyList<Vector2> TessellatedPoints,
            IReadOnlyList<double> CorrespondingKnots
        ) EquallyTessellate(NurbsCurve<Vector2> curve)
        {
            var degree = curve.Degree;
            var knots = curve.Knots;

            var uniqueKv = new List<double>();
            if (knots.Count > 0)
            {
                uniqueKv.Add(knots[0]);
                for (int k = 1; k < knots.Count; k++)
                {
                    if (!MathUtils.IsAlmostEqualTo(knots[k], uniqueKv[uniqueKv.Count - 1]))
                    {
                        uniqueKv.Add(knots[k]);
                    }
                }
            }

            var tessellatedPoints = new List<Vector2>();
            var correspondingKnots = new List<double>();

            int size = uniqueKv.Count;
            int num = 100;

            for (int i = 0; i < size - 1; i++)
            {
                double currentU = uniqueKv[i];
                double nextU = uniqueKv[i + 1];

                double step = (nextU - currentU) / (num - 1);
                for (int j = 0; j < num; j++)
                {
                    double u = currentU + step * j;
                    correspondingKnots.Add(u);
                    tessellatedPoints.Add(GetPointOnCurve(curve, u));
                }
            }

            return (tessellatedPoints, correspondingKnots);
        }

        public static double GetParamOnCurve(NurbsCurve<Vector2> curve, Vector2 givenPoint)
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            double minValue = Constants.MaxDistance;

            int maxIterations = 10;
            double paramT = Constants.DoubleEpsilon;
            double minParam = knots[0];
            double maxParam = knots[knots.Count - 1];

            var (tessellatedPoints, correspondingKnots) = EquallyTessellate(curve);

            for (int i = 0; i < tessellatedPoints.Count - 1; i++)
            {
                double currentU = correspondingKnots[i];
                double nextU = correspondingKnots[i + 1];

                Vector2 currentPoint = tessellatedPoints[i];
                Vector2 nextPoint = tessellatedPoints[i + 1];

                Vector2 diff1 = givenPoint - currentPoint;
                Vector2 vector1 =
                    diff1.LengthSquared() > MathUtil.ZeroTolerance
                        ? Vector2.Normalize(diff1)
                        : Vector2.Zero;

                Vector2 diff2 = nextPoint - currentPoint;
                Vector2 vector2 =
                    diff2.LengthSquared() > MathUtil.ZeroTolerance
                        ? Vector2.Normalize(diff2)
                        : Vector2.Zero;

                double dot = Vector2.Dot(vector1, vector2);

                Vector2 projectPoint;
                double projectU;

                if (dot < 0.0)
                {
                    projectPoint = currentPoint;
                    projectU = currentU;
                }
                else if (dot >= 1.0)
                {
                    projectPoint = nextPoint;
                    projectU = nextU;
                }
                else
                {
                    projectPoint = currentPoint + (nextPoint - currentPoint) * (float)dot;
                    projectU = currentU + (nextU - currentU) * dot;
                }

                double distance = (givenPoint - projectPoint).Length();
                if (distance < minValue)
                {
                    minValue = distance;
                    paramT = projectU;
                }
            }

            bool isClosed = IsClosed(curve);
            double a = minParam;
            double b = maxParam;

            int counters = 0;
            while (counters < maxIterations)
            {
                var derivatives = ComputeRationalCurveDerivatives(curve, 2, paramT);
                Vector2 difference = derivatives[0] - givenPoint;

                Vector2 der1 = derivatives[1];
                double f = Vector2.Dot(der1, difference);

                double condition1 = difference.Length();
                double der1Len = der1.Length();

                double denom = der1Len * condition1;
                double condition2 = MathUtils.IsZero(denom) ? 0.0 : Math.Abs(f / denom);

                if (
                    condition1 < Constants.DistanceEpsilon
                    && condition2 < Constants.DistanceEpsilon
                )
                {
                    return paramT;
                }

                Vector2 der2 = derivatives[2];
                double df = Vector2.Dot(der2, difference) + Vector2.Dot(der1, der1);

                double temp = paramT;
                if (!MathUtils.IsZero(df))
                {
                    temp = paramT - f / df;
                }

                if (!isClosed)
                {
                    if (temp < a)
                    {
                        temp = a;
                    }
                    if (temp > b)
                    {
                        temp = b;
                    }
                }
                else
                {
                    if (temp < a)
                    {
                        temp = b - (a - temp);
                    }
                    if (temp > b)
                    {
                        temp = a + (temp - b);
                    }
                }

                double condition4 = ((temp - paramT) * derivatives[1]).Length();
                if (condition4 < Constants.DistanceEpsilon)
                {
                    return paramT;
                }

                paramT = temp;
                counters++;
            }
            return paramT;
        }

        public static NurbsCurve<Vector2> Reparametrize(
            NurbsCurve<Vector2> curve,
            double alpha,
            double beta,
            double gamma,
            double delta
        )
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Argument(
                MathUtils.IsGreaterThan(alpha * delta, gamma * beta),
                "coefficient",
                "(alpha * delta - gamma * beta) must be greater than zero"
            );

            var updatedKnotVector = new double[knots.Count];
            for (int i = 0; i < knots.Count; i++)
            {
                updatedKnotVector[i] = (alpha * knots[i] + beta) / (gamma * knots[i] + delta);
            }

            var updatedControlPoints = new ControlPoint<Vector2>[controlPoints.Count];
            for (int i = 0; i < controlPoints.Count; i++)
            {
                double temp = 1.0;
                for (int j = 1; j <= degree; j++)
                {
                    double lambda = updatedKnotVector[i + j] * gamma - alpha;
                    temp = temp * lambda;
                }
                double newW = Math.Abs(controlPoints[i].Weight * temp);
                updatedControlPoints[i] = new ControlPoint<Vector2>
                {
                    Value = controlPoints[i].Value,
                    Weight = newW,
                };
            }

            return new NurbsCurve<Vector2>(degree, updatedControlPoints, updatedKnotVector);
        }

        public static NurbsCurve<Vector2> Reparametrize(
            NurbsCurve<Vector2> curve,
            double min,
            double max
        )
        {
            var knots = curve.Knots;

            if (
                MathUtils.IsAlmostEqualTo(min, knots[0])
                && MathUtils.IsAlmostEqualTo(max, knots[knots.Count - 1])
            )
            {
                return curve;
            }

            var newKnotVector = new double[knots.Count];
            double oldMin = knots[0];
            double oldMax = knots[knots.Count - 1];
            double scale = (max - min) / (oldMax - oldMin);

            for (int i = 0; i < knots.Count; i++)
            {
                newKnotVector[i] = min + (knots[i] - oldMin) * scale;
            }

            return new NurbsCurve<Vector2>(curve.Degree, curve.ControlPoints, newKnotVector);
        }

        public static bool CanComputeDerivative(NurbsCurve<Vector2> curve, double paramT)
        {
            var knots = curve.Knots;
            if (
                MathUtils.IsAlmostEqualTo(paramT, knots[0])
                || MathUtils.IsAlmostEqualTo(paramT, knots[knots.Count - 1])
            )
            {
                return true;
            }

            var pt = GetPointOnCurve(curve, paramT);
            double h = Constants.DoubleEpsilon;

            var left = (GetPointOnCurve(curve, paramT - h) - pt) / (float)-h;
            var right = (GetPointOnCurve(curve, paramT + h) - pt) / (float)h;

            return MathUtils.IsAlmostEqualTo(left, right);
        }

        public static double GetCurvature(NurbsCurve<Vector2> curve, double paramT)
        {
            var knots = curve.Knots;
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            var derivatives = ComputeRationalCurveDerivatives(curve, 2, paramT);
            var d1 = derivatives[1];
            var d2 = derivatives[2];

            double d1Length = d1.Length();

            if (MathUtils.IsAlmostEqualTo(d1Length, 1.0))
            {
                return d2.Length();
            }

            // In 2D, the cross product magnitude is abs(x1*y2 - y1*x2)
            double numerator = Math.Abs(d1.X * d2.Y - d1.Y * d2.X);
            double denominator = Math.Pow(d1Length, 3);

            if (MathUtils.IsZero(denominator))
                return 0.0;

            return numerator / denominator;
        }

        public static double GetTorsion(NurbsCurve<Vector2> curve, double paramT)
        {
            // Torsion is always zero for planar curves (2D).
            return 0.0;
        }

        public static Vector2 GetNormal(
            NurbsCurve<Vector2> curve,
            CurveNormal normalType,
            double paramT
        )
        {
            var knots = curve.Knots;
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            if (normalType == CurveNormal.Binormal)
                return Vector2.Zero;

            var derivatives = ComputeRationalCurveDerivatives(curve, 1, paramT);
            var tangent = derivatives[1];

            if (MathUtils.IsZero(tangent.LengthSquared()))
                return Vector2.Zero;

            var t = Vector2.Normalize(tangent);
            return new Vector2(-t.Y, t.X);
        }

        public static IReadOnlyList<Vector2> ProjectNormal(NurbsCurve<Vector2> curve)
        {
            // In 2D, the normal is well-defined by the tangent (rotation by 90 deg).
            // We just evaluate it at each knot.
            var knots = curve.Knots;
            int size = knots.Count;
            var normals = new Vector2[size];

            for (int i = 0; i < size; i++)
            {
                normals[i] = GetNormal(curve, CurveNormal.Normal, knots[i]);
            }

            return normals;
        }

        public static bool IsLinear(NurbsCurve<Vector2> curve)
        {
            var count = curve.ControlPoints.Count;
            if (count < 2)
                return false;
            if (count == 2)
                return true;

            float tolerance = (float)Constants.DoubleEpsilon;

            var p0 = curve.ControlPoints[0].Value;
            var v0 = Vector2.Zero;

            int i = 1;
            for (; i < count; i++)
            {
                var pi = curve.ControlPoints[i].Value;
                v0 = pi - p0;
                if (v0.LengthSquared() > tolerance)
                    break;
            }

            if (i == count)
                return true;

            v0 = Vector2.Normalize(v0);

            for (int k = 2; k < count; k++)
            {
                var pk = curve.ControlPoints[k].Value;
                var vk = pk - p0;

                // 2D Cross Product (Z component)
                float cross = v0.X * vk.Y - v0.Y * vk.X;
                if (Math.Abs(cross) > tolerance)
                    return false;
            }
            return true;
        }

        public static bool IsClosed(NurbsCurve<Vector2> curve)
        {
            var knots = curve.Knots;
            double first = knots[0];
            double end = knots[knots.Count - 1];

            var startPoint = GetPointOnCurve(curve, first);
            var endPoint = GetPointOnCurve(curve, end);

            // Geometric closure check
            if (MathUtils.IsAlmostEqualTo(startPoint, endPoint))
                return true;

            // Control point wrap-around check
            var controlPoints = curve.ControlPoints;
            int n = controlPoints.Count - 1;
            var lastCP = controlPoints[n].Value;

            int index = -1;
            for (int i = 0; i < n; i++)
            {
                var current = controlPoints[i].Value;
                if (MathUtils.IsAlmostEqualTo(lastCP, current))
                {
                    index = i;
                    break;
                }
            }

            if (index == -1)
                return false;

            if (index == 0)
                return true;

            for (int i = index; i >= 0; i--)
            {
                var current = controlPoints[i].Value;
                var another = controlPoints[n - index + i].Value;
                if (!MathUtils.IsAlmostEqualTo(current, another))
                {
                    return false;
                }
            }

            return true;
        }

        private static void ComputeDerivatives(
            NurbsCurve<Vector2> curve,
            int derivative,
            double paramT,
            Span<Vector3> result
        )
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Argument(
                derivative >= 0,
                "derivative",
                "Derivative must be greater than or equal to zero."
            );
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            int du = Math.Min(derivative, degree);
            int spanIndex = Polynomials.GetKnotSpanIndex(degree, knots, paramT);
            double[][] nders = Polynomials.BasisFunctionsDerivatives(
                spanIndex,
                degree,
                du,
                knots,
                paramT
            );

            // Fill the results span
            for (int k = 0; k <= du; k++)
            {
                result[k] = Vector3.Zero;
                for (int j = 0; j <= degree; j++)
                {
                    var cp = controlPoints[spanIndex - degree + j];
                    double w = cp.Weight;
                    // Pre-multiply weight for homogeneous coordinate
                    var pw = new Vector3(cp.Value * (float)w, (float)w);

                    result[k] += (float)nders[k][j] * pw;
                }
            }
        }
    }
}
