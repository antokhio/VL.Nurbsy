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
                    for (int j = degree, coll = 0; j > multi; j--, coll++)
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

        public static bool RemoveKnot(
            NurbsCurve<Vector3> curve,
            double removeKnot,
            int times,
            out NurbsCurve<Vector3> result
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

            for (int k = r + 1; k <= m; k++)
            {
                restKnotVector[k - times] = restKnotVector[k];
            }
            restKnotVector.RemoveRange(restKnotVector.Count - times, times);

            var updatedControlPoints = new ControlPoint<Vector3>[controlPoints.Count];
            for (int k = 0; k < controlPoints.Count; k++)
                updatedControlPoints[k] = controlPoints[k];

            var temp = new Vector4[2 * degree + 1];

            int t = 0;
            for (t = 0; t < times; t++)
            {
                int off = first - 1;
                var cpOff = updatedControlPoints[off];
                temp[0] = new Vector4(cpOff.Value * (float)cpOff.Weight, (float)cpOff.Weight);

                var cpLast = updatedControlPoints[last + 1];
                temp[last + 1 - off] = new Vector4(
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
                            new Vector4(
                                updatedControlPoints[i].Value
                                    * (float)updatedControlPoints[i].Weight,
                                (float)updatedControlPoints[i].Weight
                            )
                            - (float)(1.0 - alphai) * temp[ii - 1]
                        ) / (float)alphai;

                    temp[jj] =
                        (
                            new Vector4(
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
                            Vector4.Distance(temp[ii - 1], temp[jj + 1]),
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
                    var vecI = new Vector4(cpI.Value * (float)cpI.Weight, (float)cpI.Weight);
                    var checkVec =
                        (float)alphai * temp[ii + t + 1] + (float)(1.0 - alphai) * temp[ii - 1];
                    if (MathUtils.IsLessThanOrEqual(Vector4.Distance(vecI, checkVec), tol))
                    {
                        remflag = true;
                    }
                }

                if (!remflag)
                {
                    break;
                }

                i = first;
                j = last;

                while (j - i > t)
                {
                    var tI = temp[i - off];
                    updatedControlPoints[i] = MathUtils.IsZero(tI.W)
                        ? new ControlPoint<Vector3>(Vector3.Zero, 0)
                        : new ControlPoint<Vector3>(new Vector3(tI.X, tI.Y, tI.Z) / tI.W, tI.W);

                    var tJ = temp[j - off];
                    updatedControlPoints[j] = MathUtils.IsZero(tJ.W)
                        ? new ControlPoint<Vector3>(Vector3.Zero, 0)
                        : new ControlPoint<Vector3>(new Vector3(tJ.X, tJ.Y, tJ.Z) / tJ.W, tJ.W);

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

            int jj2 = (2 * r - s - degree) / 2;
            int ii2 = jj2;

            for (int k = 1; k < t; k++)
            {
                if (k % 2 == 1)
                    ii2++;
                else
                    jj2--;
            }

            int currJ = jj2;
            for (int k = ii2 + 1; k <= n; k++)
            {
                updatedControlPoints[currJ] = controlPoints[k];
                currJ++;
            }

            var finalCPs = new List<ControlPoint<Vector3>>(updatedControlPoints);
            finalCPs.RemoveRange(finalCPs.Count - t, t);

            result = new NurbsCurve<Vector3>(degree, finalCPs, restKnotVector);
            return true;
        }

        public static NurbsCurve<Vector3> RemoveExcessiveKnots(NurbsCurve<Vector3> curve)
        {
            var result = curve;
            var map = KnotsUtils.GetInternalKnotMultiplicityMap(curve.Knots);

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

        public static NurbsCurve<Vector3> ElevateDegree(NurbsCurve<Vector3> curve, int times)
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
            var updatedControlPoints = new List<Vector4>(
                Enumerable.Repeat(
                    new Vector4(
                        (float)Constants.MaxDistance,
                        (float)Constants.MaxDistance,
                        (float)Constants.MaxDistance,
                        1.0f
                    ),
                    moresize
                )
            );

            updatedControlPoints[0] = new Vector4(
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

            var bpts = new Vector4[degree + 1];
            for (int i = 0; i <= degree; i++)
                bpts[i] = new Vector4(
                    controlPoints[i].Value * (float)controlPoints[i].Weight,
                    (float)controlPoints[i].Weight
                );

            var nextbpts = new Vector4[degree];

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
                        alfs[k - mul - 1] = numer / (knots[a + k] - ua);
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

                var ebpts = new Vector4[degree + times + 1];
                for (int ii = lbz; ii <= ph; ii++)
                {
                    ebpts[ii] = Vector4.Zero;
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
                        bpts[j] = new Vector4(cp.Value * (float)cp.Weight, (float)cp.Weight);
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

            var finalControlPoints = new List<ControlPoint<Vector3>>();
            foreach (var h in updatedControlPoints)
            {
                if (MathUtils.IsZero(h.W))
                    finalControlPoints.Add(new ControlPoint<Vector3>(Vector3.Zero, 0));
                else
                    finalControlPoints.Add(
                        new ControlPoint<Vector3>(new Vector3(h.X, h.Y, h.Z) / h.W, h.W)
                    );
            }

            return new NurbsCurve<Vector3>(ph, finalControlPoints, updatedKnotVector);
        }

        public static bool ReduceDegree(NurbsCurve<Vector3> curve, out NurbsCurve<Vector3> result)
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
            var updatedControlPoints = new Vector4[degree];

            var cp0 = controlPoints[0];
            updatedControlPoints[0] = new Vector4(cp0.Value * (float)cp0.Weight, (float)cp0.Weight);
            var cpDeg = controlPoints[degree];
            updatedControlPoints[degree - 1] = new Vector4(
                cpDeg.Value * (float)cpDeg.Weight,
                (float)cpDeg.Weight
            );

            var homCP = new Vector4[size];
            for (int i = 0; i < size; i++)
            {
                var cp = controlPoints[i];
                homCP[i] = new Vector4(cp.Value * (float)cp.Weight, (float)cp.Weight);
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
                error = Vector4.Distance(homCP[r + 1], midP);
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
                error = Vector4.Distance(PLr, PRr);
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

            var sortedKeys = map.Keys.OrderBy(k => k).ToList();
            foreach (var u in sortedKeys)
            {
                int count = map[u] - 1;
                for (int i = 0; i < count; i++)
                    updatedKnotVector.Add(u);
            }

            var finalCPs = new List<ControlPoint<Vector3>>();
            foreach (var h in updatedControlPoints)
            {
                if (MathUtils.IsZero(h.W))
                    finalCPs.Add(new ControlPoint<Vector3>(Vector3.Zero, 0));
                else
                    finalCPs.Add(new ControlPoint<Vector3>(new Vector3(h.X, h.Y, h.Z) / h.W, h.W));
            }

            result = new NurbsCurve<Vector3>(degree - 1, finalCPs, updatedKnotVector);
            return true;
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

        private static double ComputeCurveModifyTolerance(
            IReadOnlyList<ControlPoint<Vector3>> controlPoints
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
    }
}
