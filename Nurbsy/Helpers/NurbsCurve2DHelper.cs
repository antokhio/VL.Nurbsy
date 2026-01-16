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

        public static NurbsCurve<Vector2> RefineKnotVector(
            NurbsCurve<Vector2> curve,
            IReadOnlyList<double> insertKnotElements
        )
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Argument(
                insertKnotElements != null && insertKnotElements.Count > 0,
                "insertKnotElements",
                "insertKnotElements size must be greater than zero."
            );

            int n = controlPoints.Count - 1;
            int m = n + degree + 1;
            int r = insertKnotElements.Count - 1;

            int a = Polynomials.GetKnotSpanIndex(degree, knots, insertKnotElements[0]);
            int b = Polynomials.GetKnotSpanIndex(degree, knots, insertKnotElements[r]) + 1;

            var insertedKnotVector = new double[m + r + 2];
            for (int j = 0; j <= a; j++)
            {
                insertedKnotVector[j] = knots[j];
            }
            for (int j = b + degree; j <= m; j++)
            {
                insertedKnotVector[j + r + 1] = knots[j];
            }

            var updatedControlPoints = new ControlPoint<Vector2>[n + r + 2];
            for (int j = 0; j <= a - degree; j++)
            {
                updatedControlPoints[j] = controlPoints[j];
            }
            for (int j = b - 1; j <= n; j++)
            {
                updatedControlPoints[j + r + 1] = controlPoints[j];
            }

            int i = b + degree - 1;
            int k = b + degree + r;

            for (int j = r; j >= 0; j--)
            {
                while (i > a && insertKnotElements[j] <= knots[i])
                {
                    updatedControlPoints[k - degree - 1] = controlPoints[i - degree - 1];
                    insertedKnotVector[k] = knots[i];
                    k = k - 1;
                    i = i - 1;
                }

                updatedControlPoints[k - degree - 1] = updatedControlPoints[k - degree];

                for (int l = 1; l <= degree; l++)
                {
                    int ind = k - degree + l;
                    double alpha = insertedKnotVector[k + l] - insertKnotElements[j];

                    if (MathUtils.IsAlmostEqualTo(Math.Abs(alpha), 0.0))
                    {
                        updatedControlPoints[ind - 1] = updatedControlPoints[ind];
                    }
                    else
                    {
                        alpha = alpha / (insertedKnotVector[k + l] - knots[i - degree + l]);

                        var cp1 = updatedControlPoints[ind - 1];
                        var cp2 = updatedControlPoints[ind];

                        Vector3 v1 = new Vector3(cp1.Value * (float)cp1.Weight, (float)cp1.Weight);
                        Vector3 v2 = new Vector3(cp2.Value * (float)cp2.Weight, (float)cp2.Weight);

                        Vector3 mixed = (float)alpha * v1 + (float)(1.0 - alpha) * v2;

                        if (MathUtils.IsZero(mixed.Z))
                        {
                            updatedControlPoints[ind - 1] = new ControlPoint<Vector2>(
                                Vector2.Zero,
                                0
                            );
                        }
                        else
                        {
                            updatedControlPoints[ind - 1] = new ControlPoint<Vector2>(
                                new Vector2(mixed.X, mixed.Y) / mixed.Z,
                                mixed.Z
                            );
                        }
                    }
                }
                insertedKnotVector[k] = insertKnotElements[j];
                k = k - 1;
            }

            return new NurbsCurve<Vector2>(degree, updatedControlPoints, insertedKnotVector);
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
                (alpha * delta - gamma * beta) > 0.0,
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

                double oldW = controlPoints[i].Weight;
                double newW = Math.Abs(oldW * temp);

                updatedControlPoints[i] = new ControlPoint<Vector2>(controlPoints[i].Value, newW);
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

            double currentMin = knots[0];
            double currentMax = knots[knots.Count - 1];
            double scale = (max - min) / (currentMax - currentMin);

            var newKnots = new double[knots.Count];
            for (int i = 0; i < knots.Count; i++)
            {
                newKnots[i] = min + (knots[i] - currentMin) * scale;
            }

            return new NurbsCurve<Vector2>(curve.Degree, curve.ControlPoints, newKnots);
        }

        public static NurbsCurve<Vector2> Reverse(NurbsCurve<Vector2> curve)
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            int size = knots.Count;
            var reversedKnots = new double[size];
            double min = knots[0];

            reversedKnots[0] = min;
            for (int i = 1; i < size; i++)
            {
                reversedKnots[i] = reversedKnots[i - 1] + (knots[size - i] - knots[size - i - 1]);
            }

            var reversedCPs = new List<ControlPoint<Vector2>>(controlPoints);
            reversedCPs.Reverse();

            return new NurbsCurve<Vector2>(degree, reversedCPs, reversedKnots);
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
                        temp = a;
                    if (temp > b)
                        temp = b;
                }
                else
                {
                    if (temp < a)
                        temp = b - (a - temp);
                    if (temp > b)
                        temp = a + (temp - b);
                }

                double condition4 = ((float)(temp - paramT) * derivatives[1]).Length();
                if (condition4 < Constants.DistanceEpsilon)
                {
                    return paramT;
                }

                paramT = temp;
                counters++;
            }
            return paramT;
        }

        public static bool SplitAt(
            NurbsCurve<Vector2> curve,
            double parameter,
            out NurbsCurve<Vector2> left,
            out NurbsCurve<Vector2> right
        )
        {
            int degree = curve.Degree;
            var knots = curve.Knots;

            left = default;
            right = default;

            if (parameter <= knots[degree] || parameter >= knots[knots.Count - degree - 1])
            {
                if (
                    MathUtils.IsAlmostEqualTo(parameter, knots[degree])
                    || parameter < knots[degree]
                )
                    return false;
                if (
                    MathUtils.IsAlmostEqualTo(parameter, knots[knots.Count - degree - 1])
                    || parameter > knots[knots.Count - degree - 1]
                )
                    return false;
            }

            int multi = Polynomials.GetKnotMultiplicity(knots, parameter);
            int needed = degree + 1 - multi;
            var insert = new double[needed];
            for (int k = 0; k < needed; k++)
                insert[k] = parameter;

            var tempLeft = curve;
            if (needed > 0)
            {
                tempLeft = RefineKnotVector(curve, insert);
            }

            var lKnots = tempLeft.Knots;
            var lCPs = tempLeft.ControlPoints;
            int spanIndex =
                Polynomials.GetKnotSpanIndex(tempLeft.Degree, lKnots, parameter) - degree;

            int rControlPointsCount = lCPs.Count - spanIndex;
            var rightCPs = new ControlPoint<Vector2>[rControlPointsCount];
            var rightKnots = new double[rControlPointsCount + degree + 1];

            for (int i = lCPs.Count - 1, j = rControlPointsCount - 1; j >= 0; j--, i--)
            {
                rightCPs[j] = lCPs[i];
            }

            for (int i = lKnots.Count - 1, j = rControlPointsCount + degree; j >= 0; j--, i--)
            {
                rightKnots[j] = lKnots[i];
            }

            right = new NurbsCurve<Vector2>(degree, rightCPs, rightKnots);

            var leftCPs = new ControlPoint<Vector2>[spanIndex];
            for (int i = 0; i < spanIndex; i++)
                leftCPs[i] = lCPs[i];

            var leftKnots = new double[spanIndex + degree + 1];
            for (int i = 0; i < leftKnots.Length; i++)
                leftKnots[i] = lKnots[i];

            left = new NurbsCurve<Vector2>(degree, leftCPs, leftKnots);

            return true;
        }

        public static bool Segment(
            NurbsCurve<Vector2> curve,
            double startParameter,
            double endParameter,
            out NurbsCurve<Vector2> segment
        )
        {
            segment = default;
            var degree = curve.Degree;
            var knots = curve.Knots;

            bool startAtBeginning =
                (startParameter < knots[degree])
                || MathUtils.IsAlmostEqualTo(startParameter, knots[degree]);
            bool endAtEnd =
                (endParameter > knots[knots.Count - degree - 1])
                || MathUtils.IsAlmostEqualTo(endParameter, knots[knots.Count - degree - 1]);

            if (startAtBeginning)
            {
                if (endAtEnd)
                {
                    segment = curve;
                    return true;
                }

                return SplitAt(curve, endParameter, out segment, out _);
            }
            else
            {
                if (endAtEnd)
                {
                    return SplitAt(curve, startParameter, out _, out segment);
                }
                else
                {
                    if (SplitAt(curve, startParameter, out _, out var right))
                    {
                        return SplitAt(right, endParameter, out segment, out _);
                    }
                }
            }
            return false;
        }

        public static int InsertKnot(
            NurbsCurve<Vector2> curve,
            double insertKnot,
            int times,
            out NurbsCurve<Vector2> result
        )
        {
            int degree = curve.Degree;
            var knotVector = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Argument(times > 0, nameof(times), "Times must be greater than zero.");
            Validate.Range(
                insertKnot,
                knotVector[0],
                knotVector[knotVector.Count - 1],
                nameof(insertKnot)
            );

            int knotSpanIndex = Polynomials.GetKnotSpanIndex(degree, knotVector, insertKnot);
            int originMultiplicity = Polynomials.GetKnotMultiplicity(knotVector, insertKnot);

            if (originMultiplicity + times > degree + 1)
            {
                times = degree + 1 - originMultiplicity;
            }

            if (times <= 0)
            {
                result = curve;
                return 0;
            }

            var insertedKnotVector = new List<double>(knotVector.Count + times);
            for (int i = 0; i <= knotSpanIndex; i++)
            {
                insertedKnotVector.Add(knotVector[i]);
            }
            for (int i = 1; i <= times; i++)
            {
                insertedKnotVector.Add(insertKnot);
            }
            for (int i = knotSpanIndex + 1; i < knotVector.Count; i++)
            {
                insertedKnotVector.Add(knotVector[i]);
            }

            var updatedControlPoints = new ControlPoint<Vector2>[controlPoints.Count + times];
            for (int i = 0; i <= knotSpanIndex - degree; i++)
            {
                updatedControlPoints[i] = controlPoints[i];
            }
            for (int i = knotSpanIndex - originMultiplicity; i < controlPoints.Count; i++)
            {
                updatedControlPoints[i + times] = controlPoints[i];
            }

            var temp = new Vector3[degree - originMultiplicity + 1];
            for (int i = 0; i <= degree - originMultiplicity; i++)
            {
                var cp = controlPoints[knotSpanIndex - degree + i];
                temp[i] = new Vector3(cp.Value * (float)cp.Weight, (float)cp.Weight);
            }

            int L = 0;
            for (int j = 1; j <= times; j++)
            {
                L = knotSpanIndex - degree + j;
                for (int i = 0; i <= degree - j - originMultiplicity; i++)
                {
                    double alpha =
                        (insertKnot - knotVector[L + i])
                        / (knotVector[i + knotSpanIndex + 1] - knotVector[L + i]);
                    temp[i] = (float)alpha * temp[i + 1] + (float)(1.0 - alpha) * temp[i];
                }

                var tVal = temp[0];
                updatedControlPoints[L] = MathUtils.IsZero(tVal.Z)
                    ? new ControlPoint<Vector2>(Vector2.Zero, 0)
                    : new ControlPoint<Vector2>(new Vector2(tVal.X, tVal.Y) / tVal.Z, tVal.Z);

                if (degree - j - originMultiplicity > 0)
                {
                    var tBack = temp[degree - j - originMultiplicity];
                    int idx = knotSpanIndex + times - j - originMultiplicity;
                    updatedControlPoints[idx] = MathUtils.IsZero(tBack.Z)
                        ? new ControlPoint<Vector2>(Vector2.Zero, 0)
                        : new ControlPoint<Vector2>(
                            new Vector2(tBack.X, tBack.Y) / tBack.Z,
                            tBack.Z
                        );
                }
            }

            for (int i = L + 1; i < knotSpanIndex - originMultiplicity; i++)
            {
                var tVal = temp[i - L];
                updatedControlPoints[i] = MathUtils.IsZero(tVal.Z)
                    ? new ControlPoint<Vector2>(Vector2.Zero, 0)
                    : new ControlPoint<Vector2>(new Vector2(tVal.X, tVal.Y) / tVal.Z, tVal.Z);
            }

            result = new NurbsCurve<Vector2>(
                degree,
                updatedControlPoints,
                insertedKnotVector.ToArray()
            );
            return times;
        }

        public static bool RemoveKnot(
            NurbsCurve<Vector2> curve,
            double removeKnot,
            int times,
            out NurbsCurve<Vector2> result
        )
        {
            int degree = curve.Degree;
            var knotVector = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Range(
                removeKnot,
                knotVector[0],
                knotVector[knotVector.Count - 1],
                nameof(removeKnot)
            );
            Validate.Argument(times > 0, nameof(times), "Times must be greater than zero.");

            double tol = ComputeCurveModifyTolerance(controlPoints);
            int n = controlPoints.Count - 1;

            int order = degree + 1;
            int s = Polynomials.GetKnotMultiplicity(knotVector, removeKnot);
            int r = Polynomials.GetKnotSpanIndex(degree, knotVector, removeKnot);

            int first = r - degree;
            int last = r - s;

            var restKnotVector = new List<double>(knotVector);
            int m = n + degree + 1;
            // The algorithm modifies restKnotVector. In C++: restKnotVector[k - times] = restKnotVector[k];
            // Since we need to modify in-place shifting elements, convert to array or list ops.
            // But loop runs r+1 to m. K starts after insertion point.
            // Actually it removes 'times' knots.
            for (int k = r + 1; k <= m; k++)
            {
                restKnotVector[k - times] = restKnotVector[k];
            }
            // Remove last 'times' elements
            restKnotVector.RemoveRange(restKnotVector.Count - times, times);

            // Working with homogeneous coordinates W*P
            var updatedControlPoints = new ControlPoint<Vector2>[controlPoints.Count];
            for (int k = 0; k < controlPoints.Count; k++)
                updatedControlPoints[k] = controlPoints[k];

            var temp = new Vector3[2 * degree + 1];

            int t = 0;
            for (t = 0; t < times; t++)
            {
                int off = first - 1;
                // Load temp
                var cpOff = updatedControlPoints[off];
                temp[0] = new Vector3(cpOff.Value * (float)cpOff.Weight, (float)cpOff.Weight);

                var cpLast = updatedControlPoints[last + 1];
                temp[last + 1 - off] = new Vector3(
                    cpLast.Value * (float)cpLast.Weight,
                    (float)cpLast.Weight
                );

                int i = first;
                int j = last;
                int ii = 1;
                int jj = last - off;
                bool remflag = false;

                while (j - i >= t)
                {
                    double alphai =
                        (removeKnot - knotVector[i]) / (knotVector[i + order + t] - knotVector[i]);
                    double alphaj =
                        (removeKnot - knotVector[j - t])
                        / (knotVector[j + order] - knotVector[j - t]);

                    temp[ii] =
                        (
                            new Vector3(
                                updatedControlPoints[i].Value
                                    * (float)updatedControlPoints[i].Weight,
                                (float)updatedControlPoints[i].Weight
                            )
                            - (float)(1.0 - alphai) * temp[ii - 1]
                        ) / (float)alphai;

                    temp[jj] =
                        (
                            new Vector3(
                                updatedControlPoints[j].Value
                                    * (float)updatedControlPoints[j].Weight,
                                (float)updatedControlPoints[j].Weight
                            )
                            - (float)alphaj * temp[jj + 1]
                        ) / (float)(1.0 - alphaj);

                    i++;
                    ii++;
                    j--;
                    jj--;
                }

                if (j - i < t)
                {
                    if (
                        MathUtils.IsLessThanOrEqual(
                            Vector3.Distance(temp[ii - 1], temp[jj + 1]),
                            tol
                        )
                    )
                    {
                        remflag = true;
                    }
                }
                else
                {
                    double alphai =
                        (removeKnot - knotVector[i]) / (knotVector[i + order + t] - knotVector[i]);
                    var cpI = updatedControlPoints[i];
                    var vecI = new Vector3(cpI.Value * (float)cpI.Weight, (float)cpI.Weight);
                    var checkVec =
                        (float)alphai * temp[ii + t + 1] + (float)(1.0 - alphai) * temp[ii - 1];
                    if (MathUtils.IsLessThanOrEqual(Vector3.Distance(vecI, checkVec), tol))
                    {
                        remflag = true;
                    }
                }

                if (!remflag)
                {
                    // t indicates how many were successfully removed so far (from logic flow in typical algo)
                    // but here t is loop index. If fail, we break.
                    break;
                }

                i = first;
                j = last;

                while (j - i > t)
                {
                    // Convert back from Homogeneous to ControlPoint
                    var tI = temp[i - off];
                    updatedControlPoints[i] = MathUtils.IsZero(tI.Z)
                        ? new ControlPoint<Vector2>(Vector2.Zero, 0)
                        : new ControlPoint<Vector2>(new Vector2(tI.X, tI.Y) / tI.Z, tI.Z);

                    var tJ = temp[j - off];
                    updatedControlPoints[j] = MathUtils.IsZero(tJ.Z)
                        ? new ControlPoint<Vector2>(Vector2.Zero, 0)
                        : new ControlPoint<Vector2>(new Vector2(tJ.X, tJ.Y) / tJ.Z, tJ.Z);

                    i++;
                    j--;
                }

                first--;
                last++;
            }

            if (t == 0)
            {
                result = curve;
                return false;
            }

            // Adjust control points array
            // Shift points
            int jj2 = (2 * r - s - degree) / 2;
            int ii2 = jj2;

            for (int k = 1; k < t; k++)
            {
                if (k % 2 == 1)
                    ii2++;
                else
                    jj2--;
            }

            var finalCPs = new List<ControlPoint<Vector2>>(updatedControlPoints);
            // We need to implement the final shift properly on the list or array
            // Loops k from i+1 to n. Where i is ii2.
            // updatedControlPoints[j] = controlPoints[k]
            // j starts at jj2.
            int currJ = jj2;
            for (int k = ii2 + 1; k <= n; k++)
            {
                updatedControlPoints[currJ] = controlPoints[k]; // Use ORIGINAL points for the tail?
                currJ++;
            }

            finalCPs = new List<ControlPoint<Vector2>>(updatedControlPoints);
            // Remove last t elements
            finalCPs.RemoveRange(finalCPs.Count - t, t);

            result = new NurbsCurve<Vector2>(degree, finalCPs, restKnotVector);
            return true;
        }

        public static NurbsCurve<Vector2> RemoveExcessiveKnots(NurbsCurve<Vector2> curve)
        {
            var result = curve;
            var map = KnotsUtils.GetInternalKnotMultiplicityMap(curve.Knots);

            // Iterate over map
            foreach (var kvp in map)
            {
                double u = kvp.Key;
                int count = kvp.Value;

                if (RemoveKnot(result, u, count, out var tempResult))
                {
                    result = tempResult;
                }
            }
            return result;
        }

        public static NurbsCurve<Vector2> ElevateDegree(NurbsCurve<Vector2> curve, int times)
        {
            var degree = curve.Degree;
            var knotVector = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Argument(times > 0, nameof(times), "Times must be greater than zero.");

            int n = controlPoints.Count - 1;
            int m = n + degree + 1;
            int ph = degree + times;
            int ph2 = ph / 2;

            var bezalfs = new double[degree + times + 1][];
            for (int i = 0; i < bezalfs.Length; i++)
                bezalfs[i] = new double[degree + 1];

            bezalfs[0][0] = 1.0;
            bezalfs[ph][degree] = 1.0;

            for (int i = 1; i <= ph2; i++)
            {
                double inv = 1.0 / MathUtils.Binomial(ph, i);
                int mpi = Math.Min(degree, i);

                for (int j = Math.Max(0, i - times); j <= mpi; j++)
                {
                    bezalfs[i][j] =
                        inv * MathUtils.Binomial(degree, j) * MathUtils.Binomial(times, i - j);
                }
            }

            for (int i = ph2 + 1; i <= ph - 1; i++)
            {
                int mpi = Math.Min(degree, i);
                for (int j = Math.Max(0, i - times); j <= mpi; j++)
                {
                    bezalfs[i][j] = bezalfs[ph - i][degree - j];
                }
            }

            int mh = ph;
            int kind = ph + 1;
            int r = -1;
            int a = degree;
            int b = degree + 1;
            int cind = 1;
            double ua = knotVector[0];

            int moresize = controlPoints.Count + controlPoints.Count * times;
            var updatedControlPoints = new List<Vector3>(
                Enumerable.Repeat(
                    new Vector3(
                        (float)Constants.MaxDistance,
                        (float)Constants.MaxDistance,
                        (float)Constants.MaxDistance
                    ),
                    moresize
                )
            );

            // Convert CP to Homogeneous
            updatedControlPoints[0] = new Vector3(
                controlPoints[0].Value * (float)controlPoints[0].Weight,
                (float)controlPoints[0].Weight
            );

            var updatedKnotVector = new List<double>(
                Enumerable.Repeat(Constants.MaxDistance, moresize + ph + 1)
            );
            for (int i = 0; i <= ph; i++)
            {
                updatedKnotVector[i] = ua;
            }

            var bpts = new Vector3[degree + 1];
            for (int i = 0; i <= degree; i++)
                bpts[i] = new Vector3(
                    controlPoints[i].Value * (float)controlPoints[i].Weight,
                    (float)controlPoints[i].Weight
                );

            var nextbpts = new Vector3[degree]; // degree - 1 capacity needed, but safer

            // Re-implement loop
            while (b < m)
            {
                int i = b;
                while (b < m && MathUtils.IsAlmostEqualTo(knotVector[b], knotVector[b + 1]))
                {
                    b++;
                }
                int mul = b - i + 1;
                mh += mul + times;
                double ub = knotVector[b];

                int oldr = r;
                r = degree - mul;

                int lbz = oldr > 0 ? (oldr + 2) / 2 : 1;
                int rbz = r > 0 ? ph - (r + 1) / 2 : ph;

                if (r > 0)
                {
                    double numer = ub - ua;
                    var alfs = new double[degree];
                    for (int k = degree; k > mul; k--)
                    {
                        alfs[k - mul - 1] = numer / (knotVector[a + k] - ua);
                    }
                    for (int j = 1; j <= r; j++)
                    {
                        int save = r - j;
                        int s = mul + j;
                        for (int k = degree; k >= s; k--)
                        {
                            bpts[k] =
                                (float)alfs[k - s] * bpts[k]
                                + (float)(1.0 - alfs[k - s]) * bpts[k - 1];
                        }
                        nextbpts[save] = bpts[degree];
                    }
                }

                var ebpts = new Vector3[degree + times + 1];
                for (int ii = lbz; ii <= ph; ii++)
                {
                    ebpts[ii] = Vector3.Zero;
                    int mpi = Math.Min(degree, ii);
                    for (int j = Math.Max(0, ii - times); j <= mpi; j++)
                    {
                        ebpts[ii] += (float)bezalfs[ii][j] * bpts[j];
                    }
                }

                if (oldr > 1)
                {
                    int first = kind - 2;
                    int last = kind;
                    double den = ub - ua;
                    double bet = (ub - updatedKnotVector[kind - 1]) / den;

                    for (int tr = 1; tr < oldr; tr++)
                    {
                        int ii = first;
                        int jj = last;
                        int kj = jj - kind + 1;

                        while (jj - ii > tr)
                        {
                            if (ii < cind)
                            {
                                double alf =
                                    (ub - updatedKnotVector[ii]) / (ua - updatedKnotVector[ii]);
                                updatedControlPoints[ii] =
                                    (float)alf * updatedControlPoints[ii]
                                    + (float)(1.0 - alf) * updatedControlPoints[ii - 1];
                            }

                            if (jj >= lbz)
                            {
                                if (jj - tr <= kind - ph + oldr)
                                {
                                    double gam = (ub - updatedKnotVector[jj - tr]) / den;
                                    ebpts[kj] =
                                        (float)gam * ebpts[kj] + (float)(1.0 - gam) * ebpts[kj + 1];
                                }
                                else
                                {
                                    ebpts[kj] =
                                        (float)bet * ebpts[kj] + (float)(1.0 - bet) * ebpts[kj + 1];
                                }
                            }
                            ii++;
                            jj--;
                            kj--;
                        }
                        first--;
                        last++;
                    }
                }

                if (a != degree)
                {
                    for (int ii = 0; ii < ph - oldr; ii++)
                    {
                        updatedKnotVector[kind++] = ua;
                    }
                }

                for (int j = lbz; j <= rbz; j++)
                {
                    if (cind >= updatedControlPoints.Count)
                        updatedControlPoints.Add(ebpts[j]);
                    else
                        updatedControlPoints[cind] = ebpts[j];
                    cind++;
                }

                if (b < m)
                {
                    for (int j = 0; j < r; j++)
                        bpts[j] = nextbpts[j];
                    for (int j = r; j <= degree; j++)
                    {
                        var cp = controlPoints[b - degree + j];
                        bpts[j] = new Vector3(cp.Value * (float)cp.Weight, (float)cp.Weight);
                    }
                    a = b;
                    b++;
                    ua = ub;
                }
                else
                {
                    for (int ii = 0; ii <= ph; ii++)
                        updatedKnotVector[kind + ii] = ub;
                }
            }

            // Cleanup
            for (int i = updatedControlPoints.Count - 1; i > 0; i--)
            {
                if (
                    MathUtils.IsAlmostEqualTo(
                        updatedControlPoints[i].X,
                        (float)Constants.MaxDistance
                    )
                    && MathUtils.IsAlmostEqualTo(
                        updatedControlPoints[i].Y,
                        (float)Constants.MaxDistance
                    )
                    && MathUtils.IsAlmostEqualTo(
                        updatedControlPoints[i].Z,
                        (float)Constants.MaxDistance
                    )
                )
                {
                    updatedControlPoints.RemoveAt(i);
                    continue;
                }
                break;
            }
            for (int i = updatedKnotVector.Count - 1; i > 0; i--)
            {
                if (MathUtils.IsAlmostEqualTo(updatedKnotVector[i], Constants.MaxDistance))
                {
                    updatedKnotVector.RemoveAt(i);
                    continue;
                }
                break;
            }

            var finalControlPoints = new List<ControlPoint<Vector2>>();
            foreach (var h in updatedControlPoints)
            {
                if (MathUtils.IsZero(h.Z))
                    finalControlPoints.Add(new ControlPoint<Vector2>(Vector2.Zero, 0));
                else
                    finalControlPoints.Add(
                        new ControlPoint<Vector2>(new Vector2(h.X, h.Y) / h.Z, h.Z)
                    );
            }

            return new NurbsCurve<Vector2>(ph, finalControlPoints, updatedKnotVector);
        }

        public static bool ReduceDegree(NurbsCurve<Vector2> curve, out NurbsCurve<Vector2> result)
        {
            result = default;
            int degree = curve.Degree;
            var knotVector = curve.Knots;
            var controlPoints = curve.ControlPoints;

            double tol = ComputeCurveModifyTolerance(controlPoints);
            int size = controlPoints.Count;
            bool isBezier = Validate.IsValidBezier(degree, size);

            if (!isBezier)
                return false;

            int r = (degree - 1) / 2;
            var updatedControlPoints = new Vector3[degree]; // Homogeneous

            // cp[0]
            var cp0 = controlPoints[0];
            updatedControlPoints[0] = new Vector3(cp0.Value * (float)cp0.Weight, (float)cp0.Weight);
            var cpDeg = controlPoints[degree];
            updatedControlPoints[degree - 1] = new Vector3(
                cpDeg.Value * (float)cpDeg.Weight,
                (float)cpDeg.Weight
            );

            // Need working array for homogeneous CPs of original curve
            var homCP = new Vector3[size];
            for (int i = 0; i < size; i++)
            {
                var cp = controlPoints[i];
                homCP[i] = new Vector3(cp.Value * (float)cp.Weight, (float)cp.Weight);
            }

            var alpha = new double[degree];
            for (int i = 0; i < degree; i++)
                alpha[i] = (double)i / degree;

            double error = 0.0;
            if (degree % 2 == 0)
            {
                for (int i = 1; i <= r; i++)
                {
                    updatedControlPoints[i] =
                        (homCP[i] - (float)alpha[i] * updatedControlPoints[i - 1])
                        / (float)(1.0 - alpha[i]);
                }
                for (int i = degree - 2; i > r; i--)
                {
                    updatedControlPoints[i] =
                        (homCP[i + 1] - (float)(1.0 - alpha[i + 1]) * updatedControlPoints[i + 1])
                        / (float)alpha[i + 1];
                }

                var midP = (updatedControlPoints[r] + updatedControlPoints[r + 1]) * 0.5f;
                error = Vector3.Distance(homCP[r + 1], midP);
                double c = MathUtils.Binomial(degree, r + 1);
                error = error * (c * Math.Pow(0.5, r + 1) * Math.Pow(0.5, degree - r - 1));
            }
            else
            {
                for (int i = 1; i < r; i++)
                {
                    updatedControlPoints[i] =
                        (homCP[i] - (float)alpha[i] * updatedControlPoints[i - 1])
                        / (float)(1.0 - alpha[i]);
                }
                for (int i = degree - 2; i > r; i--)
                {
                    updatedControlPoints[i] =
                        (homCP[i + 1] - (float)(1.0 - alpha[i + 1]) * updatedControlPoints[i + 1])
                        / (float)alpha[i + 1];
                }

                var PLr =
                    (homCP[r] - (float)alpha[r] * updatedControlPoints[r - 1])
                    / (float)(1.0 - alpha[r]);
                var PRr =
                    (homCP[r + 1] - (float)(1.0 - alpha[r + 1]) * updatedControlPoints[r + 1])
                    / (float)alpha[r + 1];
                updatedControlPoints[r] = (PLr + PRr) * 0.5f;
                error = Vector3.Distance(PLr, PRr);
                double maxU = (degree - Math.Sqrt(degree)) / (2 * degree);
                error =
                    error
                    * 0.5
                    * (1 - alpha[r])
                    * (
                        MathUtils.Binomial(degree, r)
                        * Math.Pow(maxU, r)
                        * Math.Pow(1 - maxU, r + 1)
                        * (1 - 2 * maxU)
                    );
            }

            if (error > tol)
                return false;

            var map = KnotsUtils.GetKnotMultiplicityMap(knotVector);
            var updatedKnotVector = new List<double>();
            // Sorted iteration is preferred
            var sortedKeys = map.Keys.OrderBy(k => k).ToList();
            foreach (var u in sortedKeys)
            {
                int count = map[u] - 1;
                for (int i = 0; i < count; i++)
                    updatedKnotVector.Add(u);
            }

            var finalCPs = new List<ControlPoint<Vector2>>();
            foreach (var h in updatedControlPoints)
            {
                if (MathUtils.IsZero(h.Z))
                    finalCPs.Add(new ControlPoint<Vector2>(Vector2.Zero, 0));
                else
                    finalCPs.Add(new ControlPoint<Vector2>(new Vector2(h.X, h.Y) / h.Z, h.Z));
            }

            result = new NurbsCurve<Vector2>(degree - 1, finalCPs, updatedKnotVector);
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

        private static double ComputeCurveModifyTolerance(
            IReadOnlyList<ControlPoint<Vector2>> controlPoints
        )
        {
            double minWeight = 1.0;
            double maxDistance = 0.0;

            for (int i = 0; i < controlPoints.Count; i++)
            {
                var cp = controlPoints[i];
                if (cp.Weight < minWeight)
                    minWeight = cp.Weight;

                // C++: temp.ToXYZ(true).Length() -> Euclidean Length of point
                double len = cp.Value.Length();
                if (len > maxDistance)
                    maxDistance = len;
            }

            return Constants.DistanceEpsilon * minWeight / (1.0 + Math.Abs(maxDistance));
        }

        public static bool Merge(
            NurbsCurve<Vector2> left,
            NurbsCurve<Vector2> right,
            out NurbsCurve<Vector2> result
        )
        {
            result = default;
            var cpL = left.ControlPoints;
            var cpR = right.ControlPoints;

            if (cpL.Count == 0 || cpR.Count == 0)
                return false;

            // Check connectivity
            if (!MathUtils.IsAlmostEqualTo(cpL[cpL.Count - 1].Value, cpR[0].Value))
            {
                return false;
            }

            int degree = Math.Max(left.Degree, right.Degree);

            // Normalize Left to [0, 1]
            var tempL = Reparametrize(left, 0.0, 1.0);
            if (degree > left.Degree)
            {
                tempL = ElevateDegree(tempL, degree - left.Degree);
                // Re-normalize after elevation (knots might change range or just structure)
                // ElevateDegree keeps range but robust to re-ensure
                tempL = Reparametrize(tempL, 0.0, 1.0);
            }

            // Normalize Right to [0, 1]
            var tempR = Reparametrize(right, 0.0, 1.0);
            if (degree > right.Degree)
            {
                tempR = ElevateDegree(tempR, degree - right.Degree);
                tempR = Reparametrize(tempR, 0.0, 1.0);
            }

            // Check if clamped (multiplicity at ends)
            // Left end
            var kL = tempL.Knots;
            int lMulti = Polynomials.GetKnotMultiplicity(kL, kL[kL.Count - 1]);
            // Right start
            var kR = tempR.Knots;
            int rMulti = Polynomials.GetKnotMultiplicity(kR, kR[0]);

            if (lMulti != degree + 1 || rMulti != degree + 1)
            {
                return false;
            }

            // Merge Control Points
            var mergedPoints = new List<ControlPoint<Vector2>>(
                tempL.ControlPoints.Count + tempR.ControlPoints.Count
            );
            mergedPoints.AddRange(tempL.ControlPoints);
            mergedPoints.AddRange(tempR.ControlPoints);

            // Merge Knots
            // L is [0, 1]. R is [0, 1] -> shift R to [1, 2]
            var shiftedR = Reparametrize(tempR, 1.0, 2.0);
            var kRShift = shiftedR.Knots;

            var mergedKnots = new List<double>(kL);
            // Append knots from R, skipping the first (degree + 1) knots which are equal to 1.0
            // Since L ends at 1.0 with (degree + 1) multiplicity, we just continue from there using R's internal + end knots.
            for (int i = degree + 1; i < kRShift.Count; i++)
            {
                mergedKnots.Add(kRShift[i]);
            }

            // Rescale back to [0, 1]
            var finalKnots = KnotsUtils.Rescale(mergedKnots, 0.0, 1.0);

            result = new NurbsCurve<Vector2>(degree, mergedPoints, finalKnots);
            return true;
        }
    }
}
