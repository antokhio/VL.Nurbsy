using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    internal static class NurbsCurve3DHelper
    {
        public static Vector3 GetPointOnCurve(NurbsCurve<Vector3> curve, double paramT)
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

            // working in Homogeneous coordinates (w*x, w*y, w*z, w)
            var temp = new Vector4[times + 1];

            for (int i = 0; i <= times; i++)
            {
                var cp = controlPoints[knotSpanIndex - degree + i];
                double w = cp.Weight;
                temp[i] = new Vector4(cp.Value * (float)w, (float)w);
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

            Vector4 result = temp[0];
            if (MathUtils.IsZero(result.W))
                return Vector3.Zero;

            return new Vector3(result.X, result.Y, result.Z) / result.W;
        }

        public static Vector3 GetPointOnCurveByCornerCut(NurbsCurve<Vector3> curve, double paramT)
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
            var temp = new Vector4[times + 1];

            for (int i = 0; i <= times; i++)
            {
                var cp = controlPoints[knotSpanIndex - degree + i];
                temp[i] = new Vector4(cp.Value * (float)cp.Weight, (float)cp.Weight);
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
            if (MathUtils.IsZero(result.W))
                return Vector3.Zero;
            return new Vector3(result.X, result.Y, result.Z) / result.W;
        }

        public static IReadOnlyList<Vector3> ComputeRationalCurveDerivatives(
            NurbsCurve<Vector3> curve,
            int derivative,
            double paramT
        )
        {
            // Homogeneous Coordinates for 3D are Vector4 (x*w, y*w, z*w, w)
            // Allocate a span-compatible array for the intermediate homogeneous derivatives
            var homogeneousDerivatives = new Vector4[derivative + 1];

            ComputeDerivatives(curve, derivative, paramT, homogeneousDerivatives);

            var derivatives = new List<Vector3>(derivative + 1);
            for (int i = 0; i <= derivative; i++)
                derivatives.Add(Vector3.Zero);

            for (int k = 0; k <= derivative; k++)
            {
                var v = new Vector3(
                    homogeneousDerivatives[k].X,
                    homogeneousDerivatives[k].Y,
                    homogeneousDerivatives[k].Z
                );

                for (int i = 1; i <= k; i++)
                {
                    double binom = MathUtils.Binomial(k, i);
                    // The weight is stored in W for 3D curve homogeneous coordinates
                    double weightDerivative = homogeneousDerivatives[i].W;

                    v -= (float)(binom * weightDerivative) * derivatives[k - i];
                }

                // Divide by the weight (0-th derivative of weight is the weight itself)
                derivatives[k] = v / homogeneousDerivatives[0].W;
            }

            return derivatives;
        }

        public static NurbsCurve<Vector3> RefineKnotVector(
            NurbsCurve<Vector3> curve,
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

            var updatedControlPoints = new ControlPoint<Vector3>[n + r + 2];
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

                        Vector4 v1 = new Vector4(cp1.Value * (float)cp1.Weight, (float)cp1.Weight);
                        Vector4 v2 = new Vector4(cp2.Value * (float)cp2.Weight, (float)cp2.Weight);

                        Vector4 mixed = (float)alpha * v1 + (float)(1.0 - alpha) * v2;

                        if (MathUtils.IsZero(mixed.W))
                        {
                            updatedControlPoints[ind - 1] = new ControlPoint<Vector3>(
                                Vector3.Zero,
                                0
                            );
                        }
                        else
                        {
                            updatedControlPoints[ind - 1] = new ControlPoint<Vector3>(
                                new Vector3(mixed.X, mixed.Y, mixed.Z) / mixed.W,
                                mixed.W
                            );
                        }
                    }
                }
                insertedKnotVector[k] = insertKnotElements[j];
                k = k - 1;
            }

            return new NurbsCurve<Vector3>(degree, updatedControlPoints, insertedKnotVector);
        }

        public static IReadOnlyList<BezierCurve<Vector3>> DecomposeToBeziers(
            NurbsCurve<Vector3> curve
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
            var beziersCP = new Vector4[bezierSize][];
            for (int i = 0; i < bezierSize; i++)
            {
                beziersCP[i] = new Vector4[degree + 1];
            }

            int a = degree;
            int b = degree + 1;
            int nb = 0;

            for (int i = 0; i <= degree; i++)
            {
                var cp = controlPoints[i];
                beziersCP[nb][i] = new Vector4(cp.Value * (float)cp.Weight, (float)cp.Weight);
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
                            beziersCP[nb][k] = new Vector4(
                                cp.Value * (float)cp.Weight,
                                (float)cp.Weight
                            );
                        }

                        a = b;
                        b += 1;
                    }
                }
            }

            var result = new List<BezierCurve<Vector3>>(bezierSize);
            for (int i = 0; i < bezierSize; i++)
            {
                var pts = new List<ControlPoint<Vector3>>(degree + 1);
                for (int k = 0; k <= degree; k++)
                {
                    var h = beziersCP[i][k];
                    double w = h.W;
                    if (MathUtils.IsZero(w))
                    {
                        pts.Add(new ControlPoint<Vector3>(Vector3.Zero, 0));
                    }
                    else
                    {
                        pts.Add(
                            new ControlPoint<Vector3>(new Vector3(h.X, h.Y, h.Z) / (float)w, w)
                        );
                    }
                }
                result.Add(new BezierCurve<Vector3>(pts));
            }

            return result;
        }

        public static (
            IReadOnlyList<Vector3> TessellatedPoints,
            IReadOnlyList<double> CorrespondingKnots
        ) EquallyTessellate(NurbsCurve<Vector3> curve)
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

            var tessellatedPoints = new List<Vector3>();
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

        public static bool CanComputeDerivative(NurbsCurve<Vector3> curve, double paramT)
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

        public static double GetCurvature(NurbsCurve<Vector3> curve, double paramT)
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

            double numerator = Vector3.Cross(d1, d2).Length();
            double denominator = Math.Pow(d1Length, 3);

            if (MathUtils.IsZero(denominator))
                return 0.0;

            return numerator / denominator;
        }

        public static double GetTorsion(NurbsCurve<Vector3> curve, double paramT)
        {
            var knots = curve.Knots;
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            // Need up to 3rd derivative
            var derivatives = ComputeRationalCurveDerivatives(curve, 3, paramT);
            var d1 = derivatives[1];
            var d2 = derivatives[2];
            var d3 = derivatives[3];

            // Formula: tau = (d1 x d2) . d3 / |d1 x d2|^2
            var cross = Vector3.Cross(d1, d2);
            double denominator = cross.LengthSquared();

            if (MathUtils.IsZero(denominator))
                return 0.0;

            double numerator = Vector3.Dot(cross, d3);
            return numerator / denominator;
        }

        public static Vector3 GetNormal(
            NurbsCurve<Vector3> curve,
            CurveNormal normalType,
            double paramT
        )
        {
            var knots = curve.Knots;
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            var derivatives = ComputeRationalCurveDerivatives(curve, 2, paramT);
            var tangent = derivatives[1];
            var der2 = derivatives[2];

            if (MathUtils.IsAlmostEqualTo(tangent.Length(), 1.0))
            {
                var lenDer2 = der2.Length();
                var curveNormal = MathUtils.IsZero(lenDer2) ? Vector3.Zero : der2 / lenDer2;

                if (normalType == CurveNormal.Normal)
                {
                    return curveNormal;
                }
                else
                {
                    return Vector3.Cross(tangent, curveNormal);
                }
            }
            else
            {
                var cross = Vector3.Cross(tangent, der2);
                var crossLen = cross.Length();

                if (MathUtils.IsZero(crossLen))
                    return Vector3.Zero;

                var b = cross / crossLen;

                if (normalType == CurveNormal.Binormal)
                {
                    return b;
                }
                else
                {
                    return Vector3.Cross(b, Vector3.Normalize(tangent));
                }
            }
        }

        public static IReadOnlyList<Vector3> ProjectNormal(NurbsCurve<Vector3> curve)
        {
            var knots = curve.Knots;
            int size = knots.Count;
            int m = size - 1;

            var Blist = new Vector3[size];

            var t0 = ComputeRationalCurveDerivatives(curve, 1, knots[0])[1];
            t0 = Vector3.Normalize(t0);

            bool flag = true;
            while (flag)
            {
                bool needReCal = false;

                // Use extension method to create initial random orthogonal
                var t0Ref = t0;
                Blist[0] = MathExtensions.RandomOrthogonal(ref t0Ref);

                for (int i = 1; i <= m; i++)
                {
                    var ti = ComputeRationalCurveDerivatives(curve, 1, knots[i])[1];
                    ti = Vector3.Normalize(ti);

                    if (MathUtils.IsZero(Vector3.Cross(ti, Blist[i - 1]).LengthSquared()))
                    {
                        needReCal = true;
                        break;
                    }

                    var dot = Vector3.Dot(Blist[i - 1], ti);
                    var bi = Blist[i - 1] - dot * ti;
                    Blist[i] = Vector3.Normalize(bi);
                }

                if (!needReCal)
                    flag = false;
            }

            if (IsClosed(curve))
            {
                Blist[m] = Blist[0];

                var Baver = new Vector3[size];
                Baver[m] = Blist[0];

                for (int i = m - 1; i >= 1; i--)
                {
                    var ti = ComputeRationalCurveDerivatives(curve, 1, knots[i])[1];
                    ti = Vector3.Normalize(ti);

                    var prev = Baver[i + 1];
                    var dot = Vector3.Dot(prev, ti);
                    var bi = prev - dot * ti;
                    Baver[i] = Vector3.Normalize(bi);
                }

                for (int i = 1; i < m; i++)
                {
                    var avg = (Blist[i] + Baver[i]) * 0.5f;
                    Blist[i] = Vector3.Normalize(avg);
                }
            }

            return Blist;
        }

        public static bool IsLinear(NurbsCurve<Vector3> curve)
        {
            var count = curve.ControlPoints.Count;
            if (count < 2)
                return false;
            if (count == 2)
                return true;

            float tolerance = (float)Constants.DoubleEpsilon;

            var p0 = curve.ControlPoints[0].Value;
            var v0 = Vector3.Zero;

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

            v0 = Vector3.Normalize(v0);

            for (int k = 2; k < count; k++)
            {
                var pk = curve.ControlPoints[k].Value;
                var vk = pk - p0;

                var cross = Vector3.Cross(v0, vk);
                if (cross.LengthSquared() > tolerance * tolerance)
                    return false;
            }
            return true;
        }

        public static bool IsClosed(NurbsCurve<Vector3> curve)
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

        public static NurbsCurve<Vector3> Reparametrize(
            NurbsCurve<Vector3> curve,
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

            var updatedControlPoints = new ControlPoint<Vector3>[controlPoints.Count];
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

                updatedControlPoints[i] = new ControlPoint<Vector3>(controlPoints[i].Value, newW);
            }

            return new NurbsCurve<Vector3>(degree, updatedControlPoints, updatedKnotVector);
        }

        public static NurbsCurve<Vector3> Reparametrize(
            NurbsCurve<Vector3> curve,
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

            return new NurbsCurve<Vector3>(curve.Degree, curve.ControlPoints, newKnots);
        }

        public static NurbsCurve<Vector3> Reverse(NurbsCurve<Vector3> curve)
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

            var reversedCPs = new List<ControlPoint<Vector3>>(controlPoints);
            reversedCPs.Reverse();

            return new NurbsCurve<Vector3>(degree, reversedCPs, reversedKnots);
        }

        public static double GetParamOnCurve(NurbsCurve<Vector3> curve, Vector3 givenPoint)
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

                Vector3 currentPoint = tessellatedPoints[i];
                Vector3 nextPoint = tessellatedPoints[i + 1];

                Vector3 diff1 = givenPoint - currentPoint;
                Vector3 vector1 =
                    diff1.LengthSquared() > MathUtil.ZeroTolerance
                        ? Vector3.Normalize(diff1)
                        : Vector3.Zero;

                Vector3 diff2 = nextPoint - currentPoint;
                Vector3 vector2 =
                    diff2.LengthSquared() > MathUtil.ZeroTolerance
                        ? Vector3.Normalize(diff2)
                        : Vector3.Zero;

                double dot = Vector3.Dot(vector1, vector2);

                Vector3 projectPoint;
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
                Vector3 difference = derivatives[0] - givenPoint;
                Vector3 der1 = derivatives[1];
                double f = Vector3.Dot(der1, difference);

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

                Vector3 der2 = derivatives[2];
                double df = Vector3.Dot(der2, difference) + Vector3.Dot(der1, der1);

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

        public static bool SplitAt(
            NurbsCurve<Vector3> curve,
            double parameter,
            out NurbsCurve<Vector3> left,
            out NurbsCurve<Vector3> right
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
            var rightCPs = new ControlPoint<Vector3>[rControlPointsCount];
            var rightKnots = new double[rControlPointsCount + degree + 1];

            for (int i = lCPs.Count - 1, j = rControlPointsCount - 1; j >= 0; j--, i--)
            {
                rightCPs[j] = lCPs[i];
            }

            for (int i = lKnots.Count - 1, j = rControlPointsCount + degree; j >= 0; j--, i--)
            {
                rightKnots[j] = lKnots[i];
            }

            right = new NurbsCurve<Vector3>(degree, rightCPs, rightKnots);

            var leftCPs = new ControlPoint<Vector3>[spanIndex];
            for (int i = 0; i < spanIndex; i++)
                leftCPs[i] = lCPs[i];

            var leftKnots = new double[spanIndex + degree + 1];
            for (int i = 0; i < leftKnots.Length; i++)
                leftKnots[i] = lKnots[i];

            left = new NurbsCurve<Vector3>(degree, leftCPs, leftKnots);

            return true;
        }

        public static bool Segment(
            NurbsCurve<Vector3> curve,
            double startParameter,
            double endParameter,
            out NurbsCurve<Vector3> segment
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

        private static void ComputeDerivatives(
            NurbsCurve<Vector3> curve,
            int derivative,
            double paramT,
            Span<Vector4> result
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
                result[k] = Vector4.Zero;
                for (int j = 0; j <= degree; j++)
                {
                    var cp = controlPoints[spanIndex - degree + j];
                    double w = cp.Weight;
                    // Pre-multiply weight for homogeneous coordinate
                    var pw = new Vector4(cp.Value * (float)w, (float)w);

                    result[k] += (float)nders[k][j] * pw;
                }
            }
        }
    }
}
