using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    internal static class NurbsCurve3DHelper
    {
        public static NurbsCurve<Vector3> CreateLine(Vector3 start, Vector3 end)
        {
            Validate.Argument(
                !MathUtils.IsAlmostEqualTo(start, end),
                nameof(end),
                "start must not be equal to end."
            );

            int degree = 1;
            var controlPoints = new ControlPoint<Vector3>[]
            {
                new ControlPoint<Vector3>(start, 1.0),
                new ControlPoint<Vector3>(end, 1.0),
            };
            var knots = new double[] { 0.0, 0.0, 1.0, 1.0 };

            return new NurbsCurve<Vector3>(degree, controlPoints, knots);
        }

        public static NurbsCurve<Vector3> CreateCubicHermite(
            IReadOnlyList<Vector3> throughPoints,
            IReadOnlyList<Vector3> tangents
        )
        {
            int n = throughPoints.Count;
            Validate.Argument(
                n > 3,
                nameof(throughPoints),
                "ThroughPoints size must be greater than three."
            );
            Validate.Argument(
                n == tangents.Count,
                nameof(tangents),
                "Tangents size must be equal to throughPoints size."
            );

            var startPoint = throughPoints[0];
            var endPoint = throughPoints[n - 1];

            var uk = Interpolation.GetChordParameterization(throughPoints);

            bool isCyclePoint = MathUtils.IsAlmostEqualTo(startPoint, endPoint);
            var startTangent = tangents[0];
            var endTangent = tangents[n - 1];
            bool isCycleTangent = MathUtils.IsAlmostEqualTo(startTangent, endTangent);

            int kn = 2 * n;
            int knotSize = kn + 4;
            var knotVector = new double[knotSize];

            for (int i = 2, j = 0; i < kn + 2; i += 2, j++)
            {
                knotVector[i] = uk[j];
                knotVector[i + 1] = uk[j];
            }

            if (isCyclePoint && isCycleTangent)
            {
                knotVector[0] = knotVector[1] = uk[0] - (uk[n - 1] - uk[n - 2]);
                knotVector[kn + 2] = knotVector[kn + 3] = uk[n - 1] + uk[1] - uk[0];
            }
            else if (isCyclePoint && !isCycleTangent)
            {
                knotVector[0] = uk[0] - (uk[n - 1] - uk[n - 2]);
                knotVector[1] = knotVector[2];
                knotVector[kn + 2] = knotVector[kn];
                knotVector[kn + 3] = uk[n - 1] + uk[1] - uk[0];
            }
            else
            {
                knotVector[0] = knotVector[1] = knotVector[2];
                knotVector[kn + 2] = knotVector[kn + 3] = knotVector[kn];
            }

            var controlPoints = new ControlPoint<Vector3>[kn];
            for (int j = 0, coef = 0; j < kn; j += 2, coef++)
            {
                double i1 = knotVector[j + 3] - knotVector[j + 1];
                double i2 = knotVector[j + 4] - knotVector[j + 2];

                controlPoints[j] = new ControlPoint<Vector3>(
                    throughPoints[coef] - (float)(i1 / 3.0) * tangents[coef],
                    1.0
                );
                controlPoints[j + 1] = new ControlPoint<Vector3>(
                    throughPoints[coef] + (float)(i2 / 3.0) * tangents[coef],
                    1.0
                );
            }

            return new NurbsCurve<Vector3>(degree: 3, controlPoints, knotVector);
        }

        public static NurbsCurve<Vector3> CreateArc(
            Vector3 center,
            Vector3 xAxis,
            Vector3 yAxis,
            double startRad,
            double endRad,
            double xRadius,
            double yRadius
        )
        {
            Validate.Argument(
                !MathUtils.IsZero(xAxis, Constants.DoubleEpsilon),
                nameof(xAxis),
                "xAxis must not be zero vector."
            );
            Validate.Argument(
                !MathUtils.IsZero(yAxis, Constants.DoubleEpsilon),
                nameof(yAxis),
                "yAxis must not be zero vector."
            );
            Validate.Argument(
                MathUtils.IsGreaterThan(endRad, startRad),
                nameof(endRad),
                "endRad must be greater than startRad."
            );

            double theta = endRad - startRad;
            Validate.Range(theta, 0, 2 * Constants.Pi, nameof(theta));
            Validate.Argument(
                MathUtils.IsGreaterThan(xRadius, 0.0),
                nameof(xRadius),
                "xRadius must be greater than zero."
            );
            Validate.Argument(
                MathUtils.IsGreaterThan(yRadius, 0.0),
                nameof(yRadius),
                "yRadius must be greater than zero."
            );

            int narcs = 0;
            if (MathUtils.IsLessThanOrEqual(theta, Constants.Pi / 2.0))
            {
                narcs = 1;
            }
            else
            {
                if (MathUtils.IsLessThanOrEqual(theta, Constants.Pi))
                {
                    narcs = 2;
                }
                else if (MathUtils.IsLessThanOrEqual(theta, 3 * Constants.Pi / 2.0))
                {
                    narcs = 3;
                }
                else
                {
                    narcs = 4;
                }
            }

            double dtheta = theta / narcs;
            int n = 2 * narcs;
            int degree = 2;

            var controlPoints = new ControlPoint<Vector3>[n + 1];
            var knotVector = new double[n + degree + 2];

            double w1 = Math.Cos(dtheta / 2.0);
            Vector3 nX = Vector3.Normalize(xAxis);
            Vector3 nY = Vector3.Normalize(yAxis);

            Vector3 P0 =
                center
                + (float)(xRadius * Math.Cos(startRad)) * nX
                + (float)(yRadius * Math.Sin(startRad)) * nY;
            Vector3 T0 = (float)-Math.Sin(startRad) * nX + (float)Math.Cos(startRad) * nY;

            controlPoints[0] = new ControlPoint<Vector3>(P0, 1.0);

            int index = 0;
            double angle = startRad;

            for (int i = 1; i <= narcs; i++)
            {
                angle += dtheta;
                Vector3 P2 =
                    center
                    + (float)(xRadius * Math.Cos(angle)) * nX
                    + (float)(yRadius * Math.Sin(angle)) * nY;
                controlPoints[index + 2] = new ControlPoint<Vector3>(P2, 1.0);

                Vector3 T2 = (float)-Math.Sin(angle) * nX + (float)Math.Cos(angle) * nY;

                var type = Intersection.ComputeRays(
                    P0,
                    T0,
                    P2,
                    T2,
                    out double param0,
                    out double param2,
                    out Vector3 P1
                );

                if (type != CurveCurveIntersectionType.Intersecting)
                {
                    throw new InvalidOperationException(
                        "Failed to compute arc control points (tangents do not intersect)."
                    );
                }

                controlPoints[index + 1] = new ControlPoint<Vector3>(P1, w1);
                index = index + 2;

                if (i < narcs)
                {
                    P0 = P2;
                    T0 = T2;
                }
            }

            int j = 2 * narcs + 1;
            for (int i = 0; i < 3; i++)
            {
                knotVector[i] = 0.0;
                knotVector[i + j] = 1.0;
            }

            switch (narcs)
            {
                case 1:
                    break;
                case 2:
                    knotVector[3] = knotVector[4] = 0.5;
                    break;
                case 3:
                    knotVector[3] = knotVector[4] = 1.0 / 3.0;
                    knotVector[5] = knotVector[6] = 2.0 / 3.0;
                    break;
                case 4:
                    knotVector[3] = knotVector[4] = 0.25;
                    knotVector[5] = knotVector[6] = 0.5;
                    knotVector[7] = knotVector[8] = 0.75;
                    break;
            }

            return new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
        }

        public static bool CreateOneConicArc(
            Vector3 start,
            Vector3 startTangent,
            Vector3 end,
            Vector3 endTangent,
            Vector3 pointOnConic,
            out Vector3 projectPoint,
            out double projectPointWeight
        )
        {
            Validate.Argument(
                !MathUtils.IsZero(startTangent),
                nameof(startTangent),
                "StartTangent must not be zero vector."
            );
            Validate.Argument(
                !MathUtils.IsZero(endTangent),
                nameof(endTangent),
                "EndTangent must not be zero vector."
            );

            projectPoint = default;
            projectPointWeight = 0.0;

            var type = Intersection.ComputeRays(
                start,
                startTangent,
                end,
                endTangent,
                out _,
                out _,
                out Vector3 point
            );

            Vector3 pDiff = end - start;

            if (type == CurveCurveIntersectionType.Intersecting)
            {
                Vector3 v1p = pointOnConic - point;
                type = Intersection.ComputeRays(
                    point,
                    v1p,
                    start,
                    pDiff,
                    out double alf0,
                    out double alf2,
                    out _
                );

                if (type == CurveCurveIntersectionType.Intersecting)
                {
                    double a = Math.Sqrt(alf2 / (1.0 - alf2));
                    double u = a / (1.0 + a);

                    double dot1 = Vector3.Dot(pointOnConic - start, point - pointOnConic);
                    double dot2 = Vector3.Dot(pointOnConic - end, point - pointOnConic);
                    double dotDen = Vector3.Dot(point - pointOnConic, point - pointOnConic);

                    double num = (1.0 - u) * (1.0 - u) * dot1 + u * u * dot2;
                    double den = 2.0 * u * (1.0 - u) * dotDen;

                    projectPoint = point;
                    projectPointWeight = num / den;
                    return true;
                }
            }
            else if (type == CurveCurveIntersectionType.Parallel)
            {
                type = Intersection.ComputeRays(
                    pointOnConic,
                    startTangent,
                    start,
                    pDiff,
                    out double alf0,
                    out double alf2,
                    out _
                );

                if (type == CurveCurveIntersectionType.Intersecting)
                {
                    double a = Math.Sqrt(alf2 / (1.0 - alf2));
                    double u = a / (1.0 + a);
                    double b = 2.0 * u * (1.0 - u);
                    b = -alf0 * (1.0 - b) / b;

                    projectPoint = startTangent * (float)b;
                    projectPointWeight = 0.0;
                    return true;
                }
            }
            return false;
        }

        public static bool CreateOpenConic(
            Vector3 start,
            Vector3 startTangent,
            Vector3 end,
            Vector3 endTangent,
            Vector3 pointOnConic,
            out NurbsCurve<Vector3> curve
        )
        {
            curve = default;
            Validate.Argument(
                !MathUtils.IsZero(startTangent),
                nameof(startTangent),
                "StartTangent must not be zero vector."
            );
            Validate.Argument(
                !MathUtils.IsZero(endTangent),
                nameof(endTangent),
                "EndTangent must not be zero vector."
            );

            if (
                !CreateOneConicArc(
                    start,
                    startTangent,
                    end,
                    endTangent,
                    pointOnConic,
                    out var P1,
                    out var w1
                )
            )
            {
                return false;
            }

            int nsegs = 0;
            if (MathUtils.IsLessThanOrEqual(w1, -1.0))
            {
                return false;
            }

            if (MathUtils.IsGreaterThanOrEqual(w1, 1.0))
            {
                nsegs = 1;
            }
            else
            {
                var v1 = Vector3.Normalize(P1 - start);
                var v2 = Vector3.Normalize(end - P1);

                // 3D Angle via DotProduct and Acos
                double dot = Vector3.Dot(v1, v2);
                if (dot > 1.0)
                    dot = 1.0;
                if (dot < -1.0)
                    dot = -1.0;
                double rad = Math.Acos(dot);

                if (MathUtils.IsGreaterThan(w1, 0.0) && rad > MathUtil.DegreesToRadians(60))
                {
                    nsegs = 1;
                }
                else if (MathUtils.IsLessThan(w1, 0.0) && rad > MathUtil.DegreesToRadians(90))
                {
                    nsegs = 4;
                }
                else
                {
                    nsegs = 2;
                }
            }

            int n = 2 * nsegs;
            int j = 2 * nsegs + 1;

            int degree = 2;
            var knotVector = new double[j + degree + 1];
            var controlPoints = new ControlPoint<Vector3>[n + 1];

            for (int i = 0; i < 3; i++)
            {
                knotVector[i] = 0.0;
                knotVector[i + j] = 1.0;
            }

            controlPoints[0] = new ControlPoint<Vector3>(start, 1.0);
            controlPoints[n] = new ControlPoint<Vector3>(end, 1.0);

            if (nsegs == 1)
            {
                controlPoints[1] = new ControlPoint<Vector3>(P1, w1);
                curve = new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
                return true;
            }

            SplitArc(start, P1, w1, end, out var Q1, out var S, out var R1, out var wqr);

            if (nsegs == 2)
            {
                controlPoints[2] = new ControlPoint<Vector3>(S, 1.0);
                controlPoints[1] = new ControlPoint<Vector3>(Q1, wqr);
                controlPoints[3] = new ControlPoint<Vector3>(R1, wqr);

                knotVector[3] = knotVector[4] = 0.5;
                curve = new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
                return true;
            }

            if (nsegs == 4)
            {
                controlPoints[4] = new ControlPoint<Vector3>(S, 1.0);
                w1 = wqr;

                SplitArc(start, Q1, w1, S, out var HQ1, out var HS, out var HR1, out wqr);
                controlPoints[2] = new ControlPoint<Vector3>(HS, 1.0);
                controlPoints[1] = new ControlPoint<Vector3>(HQ1, wqr);
                controlPoints[3] = new ControlPoint<Vector3>(HR1, wqr);

                SplitArc(S, R1, w1, end, out HQ1, out HS, out HR1, out wqr);
                controlPoints[6] = new ControlPoint<Vector3>(HS, 1.0);
                controlPoints[5] = new ControlPoint<Vector3>(HQ1, wqr);
                controlPoints[7] = new ControlPoint<Vector3>(HR1, wqr);

                for (int i = 0; i < 2; i++)
                {
                    knotVector[i + 3] = 0.25;
                    knotVector[i + 5] = 0.5;
                    knotVector[i + 7] = 0.75;
                }
                curve = new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
                return true;
            }
            return false;
        }

        public static NurbsCurve<Vector3> GlobalInterpolation(
            int degree,
            IReadOnlyList<Vector3> throughPoints,
            IReadOnlyList<double> parameters = null
        )
        {
            Validate.Argument(
                degree >= 0 && degree <= Constants.NURBSMaxDegree,
                nameof(degree),
                "Degree must be greater than or equal zero and not exceed the maximun degree."
            );
            Validate.Argument(
                throughPoints.Count > degree,
                nameof(throughPoints),
                "ThroughPoints size must be greater than degree."
            );

            int size = throughPoints.Count;
            int n = size - 1;

            IReadOnlyList<double> uk;
            if (parameters == null || parameters.Count == 0)
            {
                uk = Interpolation.GetChordParameterization(throughPoints);
            }
            else
            {
                Validate.Argument(
                    parameters.Count == size,
                    nameof(parameters),
                    "Params size must be equal to throughPoints size."
                );
                uk = parameters;
            }

            var knotVector = Interpolation.AverageKnotVector(degree, uk);

            var A = new double[size][];
            for (int i = 0; i < size; i++)
                A[i] = new double[size];

            for (int i = 1; i < n; i++)
            {
                int spanIndex = Polynomials.GetKnotSpanIndex(degree, knotVector, uk[i]);
                var basis = Polynomials.BasisFunctions(spanIndex, degree, knotVector, uk[i]);

                for (int j = 0; j <= degree; j++)
                {
                    A[i][spanIndex - degree + j] = basis[j];
                }
            }
            A[0][0] = 1.0;
            A[n][n] = 1.0;

            var right = new double[size][];
            for (int i = 0; i < size; i++)
            {
                right[i] = new double[]
                {
                    throughPoints[i].X,
                    throughPoints[i].Y,
                    throughPoints[i].Z,
                };
            }

            var result = MathUtils.SolveLinearSystem(A, right);

            var controlPoints = new ControlPoint<Vector3>[size];
            for (int i = 0; i < size; i++)
            {
                controlPoints[i] = new ControlPoint<Vector3>(
                    new Vector3((float)result[i][0], (float)result[i][1], (float)result[i][2]),
                    1.0
                );
            }

            return new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
        }

        public static NurbsCurve<Vector3> GlobalInterpolation(
            int degree,
            IReadOnlyList<Vector3> throughPoints,
            IReadOnlyList<Vector3> tangents,
            double tangentFactor = 1.0
        )
        {
            Validate.Argument(degree > 0, nameof(degree), "Degree must be greater than zero.");
            Validate.Argument(
                throughPoints.Count > degree,
                nameof(throughPoints),
                "ThroughPoints size must be greater than degree."
            );
            Validate.Argument(
                tangentFactor > 0.0,
                nameof(tangentFactor),
                "TangentFactor must be greater than zero."
            );

            var unitTangents = new Vector3[tangents.Count];
            for (int i = 0; i < tangents.Count; i++)
            {
                unitTangents[i] = Vector3.Normalize(tangents[i]);
            }

            int size = throughPoints.Count;
            int n = 2 * size;

            var controlPoints = new ControlPoint<Vector3>[n];
            var knotVector = new double[n + degree + 1];

            double d = Interpolation.GetTotalChordLength(throughPoints);
            var uk = Interpolation.GetChordParameterization(throughPoints);

            switch (degree)
            {
                case 2:
                {
                    for (int i = 0; i <= degree; i++)
                    {
                        knotVector[i] = 0.0;
                        knotVector[knotVector.Length - 1 - i] = 1.0;
                    }
                    for (int i = 0; i < size - 1; i++)
                    {
                        knotVector[2 * i + degree] = uk[i];
                        knotVector[2 * i + degree + 1] = (uk[i] + uk[i + 1]) / 2.0;
                    }
                    break;
                }
                case 3:
                {
                    for (int i = 0; i <= degree; i++)
                    {
                        knotVector[i] = 0.0;
                        knotVector[knotVector.Length - 1 - i] = 1.0;
                    }
                    for (int i = 1; i < size - 1; i++)
                    {
                        knotVector[degree + 2 * i] = (2 * uk[i] + uk[i + 1]) / 3.0;
                        knotVector[degree + 2 * i + 1] = (uk[i] + 2 * uk[i + 1]) / 3.0;
                    }
                    knotVector[4] = uk[1] / 2.0;
                    knotVector[knotVector.Length - degree - 2] = (uk[size - 1] + 1.0) / 2.0;
                    break;
                }
                default:
                {
                    var uk2 = new double[2 * size];
                    for (int i = 0; i < size - 1; i++)
                    {
                        uk2[2 * i] = uk[i];
                        uk2[2 * i + 1] = (uk[i] + uk[i + 1]) / 2.0;
                    }
                    uk2[uk2.Length - 2] = (uk2[uk2.Length - 1] + uk2[uk2.Length - 3]) / 2.0;
                    knotVector = Interpolation.AverageKnotVector(degree, uk2);
                    break;
                }
            }

            var A = new double[n][];
            for (int i = 0; i < n; i++)
                A[i] = new double[n];

            for (int i = 1; i < size - 1; i++)
            {
                int spanIndex = Polynomials.GetKnotSpanIndex(degree, knotVector, uk[i]);
                var basis = Polynomials.BasisFunctions(spanIndex, degree, knotVector, uk[i]);
                var derBasis = Polynomials.BasisFunctionsDerivatives(
                    spanIndex,
                    degree,
                    1,
                    knotVector,
                    uk[i]
                );

                for (int j = 0; j <= degree; j++)
                {
                    A[2 * i][spanIndex - degree + j] = basis[j];
                    A[2 * i + 1][spanIndex - degree + j] = derBasis[1][j];
                }
            }

            A[0][0] = 1.0;
            A[1][0] = -1.0;
            A[1][1] = 1.0;
            A[n - 2][n - 2] = -1.0;
            A[n - 2][n - 1] = 1.0;
            A[n - 1][n - 1] = 1.0;

            var right = new double[n][];
            for (int i = 0; i < size; i++)
            {
                right[2 * i] = new double[]
                {
                    throughPoints[i].X,
                    throughPoints[i].Y,
                    throughPoints[i].Z,
                };
                right[2 * i + 1] = new double[]
                {
                    unitTangents[i].X * d,
                    unitTangents[i].Y * d,
                    unitTangents[i].Z * d,
                };
            }

            double d0 = knotVector[degree + 1] / degree;
            double dn = (1 - knotVector[knotVector.Length - degree - 2]) / degree;

            Vector3 dp0 = unitTangents[0];
            Vector3 dpn = unitTangents[size - 1];
            Vector3 qpn = throughPoints[size - 1];

            right[1][0] = d0 * dp0.X * d;
            right[1][1] = d0 * dp0.Y * d;
            right[1][2] = d0 * dp0.Z * d;

            right[n - 2][0] = dn * dpn.X * d;
            right[n - 2][1] = dn * dpn.Y * d;
            right[n - 2][2] = dn * dpn.Z * d;

            right[n - 1][0] = qpn.X;
            right[n - 1][1] = qpn.Y;
            right[n - 1][2] = qpn.Z;

            var result = MathUtils.SolveLinearSystem(A, right);

            for (int i = 0; i < result.Length; i++)
            {
                controlPoints[i] = new ControlPoint<Vector3>(
                    new Vector3((float)result[i][0], (float)result[i][1], (float)result[i][2]),
                    1.0
                );
            }

            return new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
        }

        public static bool CubicLocalInterpolation(
            IReadOnlyList<Vector3> throughPoints,
            out NurbsCurve<Vector3> curve
        )
        {
            curve = default;
            Validate.Argument(
                throughPoints.Count > 0,
                nameof(throughPoints),
                "ThroughPoints size must be greater than zero."
            );

            if (!Interpolation.TryComputeTangents(throughPoints, out var tangents))
            {
                return false;
            }

            int size = throughPoints.Count;
            int n = size - 1;
            int degree = 3;

            var uk = new double[size];
            uk[0] = 0.0;

            var tempControlPoints = new List<ControlPoint<Vector3>>(2 * n);

            for (int k = 0; k < n; k++)
            {
                Vector3 t0 = tangents[k];
                Vector3 t3 = tangents[k + 1];
                Vector3 p0 = throughPoints[k];
                Vector3 p3 = throughPoints[k + 1];

                double a = 16.0 - (t0 + t3).LengthSquared();
                double b = 12.0 * Vector3.Dot(p3 - p0, t0 + t3);
                double c = -36.0 * (p3 - p0).LengthSquared();

                double det = b * b - 4.0 * a * c;
                if (det < 0)
                    det = 0;
                double alpha = (-b + Math.Sqrt(det)) / (2.0 * a);

                var pk1 = p0 + (float)(alpha / 3.0) * t0;
                var pk2 = p3 - (float)(alpha / 3.0) * t3;

                uk[k + 1] = uk[k] + 3.0 * (pk1 - p0).Length();

                tempControlPoints.Add(new ControlPoint<Vector3>(pk1, 1.0));
                tempControlPoints.Add(new ControlPoint<Vector3>(pk2, 1.0));
            }

            int kvSize = 2 * degree + 2 * n;
            var knotVector = new double[kvSize];

            for (int i = 0; i <= degree; i++)
            {
                knotVector[i] = 0.0;
            }

            double maxUk = uk[n];
            if (MathUtils.IsZero(maxUk))
                maxUk = 1.0;

            int currentKnot = degree + 1;
            for (int i = 1; i < n; i++)
            {
                double val = uk[i] / maxUk;
                knotVector[currentKnot++] = val;
                knotVector[currentKnot++] = val;
            }

            for (int i = currentKnot; i < kvSize; i++)
            {
                knotVector[i] = 1.0;
            }

            var controlPoints = new ControlPoint<Vector3>[2 * n + 2];
            controlPoints[0] = new ControlPoint<Vector3>(throughPoints[0], 1.0);
            for (int i = 0; i < tempControlPoints.Count; i++)
            {
                controlPoints[i + 1] = tempControlPoints[i];
            }
            controlPoints[2 * n + 1] = new ControlPoint<Vector3>(throughPoints[n], 1.0);

            curve = new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
            return true;
        }

        public static bool LeastSquaresApproximation(
            int degree,
            IReadOnlyList<Vector3> throughPoints,
            int controlPointsCount,
            out NurbsCurve<Vector3> curve
        )
        {
            curve = default;
            Validate.Argument(
                degree >= 0 && degree <= Constants.NURBSMaxDegree,
                nameof(degree),
                "Degree must be greater than or equal zero and not exceed the maximum degree."
            );
            Validate.Argument(
                controlPointsCount > degree,
                nameof(controlPointsCount),
                "ControlPointsCount must be greater than degree."
            );

            int n = controlPointsCount - 1;
            int m = throughPoints.Count - 1;

            if (m < degree)
                return false;

            var uk = Interpolation.GetChordParameterization(throughPoints);
            var knotVector = Interpolation.ComputeKnotVector(degree, controlPointsCount, uk);

            var start = new int[n - 1];
            var end = new int[n - 1];
            var index = new int[m - 1];

            var B = new double[m - 1][];
            for (int i = 0; i < m - 1; i++)
                B[i] = new double[degree + 1];

            int dim = n - 1;
            if (dim <= 0)
            {
                if (controlPointsCount == 2)
                {
                    curve = CreateLine(throughPoints[0], throughPoints[m]);
                    return true;
                }
                return false;
            }

            var NTN = new double[dim][];
            for (int i = 0; i < dim; i++)
                NTN[i] = new double[dim];

            for (int i = 0; i <= Math.Min(degree - 1, dim - 1); i++)
            {
                start[i] = 0;
            }
            end[0] = -2;

            int rj = degree;
            int sj = degree - 1;
            int ej = -2;

            for (int i = 1; i <= m - 1; i++)
            {
                int j = Polynomials.GetKnotSpanIndex(degree, knotVector, uk[i]);
                var N = Polynomials.BasisFunctions(j, degree, knotVector, uk[i]);

                int l = (j == degree) ? 1 : 0;
                int hk = (j == degree || j == n) ? degree - 1 : degree;

                for (int k = 0; k <= hk; k++)
                {
                    B[i - 1][k] = N[l + k];
                }

                index[i - 1] = Math.Max(0, j - degree - 1);

                if (j > rj)
                {
                    for (int k = 1; k <= j - rj; k++)
                    {
                        sj++;
                        ej++;
                        if (sj < dim)
                            start[sj] = i - 1;
                        if (ej >= 0)
                            end[ej] = i - 2;
                    }
                    rj = j;
                }
            }

            if (sj < dim - 1 || end[0] == -1)
                return false;

            for (int i = Math.Max(0, ej + 1); i < dim; i++)
            {
                end[i] = m - 2;
            }

            for (int i = 0; i < dim; i++)
            {
                int lj = Math.Max(0, i - degree);
                int hj = Math.Min(dim - 1, i + degree);
                for (int j = lj; j <= hj; j++)
                {
                    int lk = Math.Max(start[i], start[j]);
                    int hk = Math.Min(end[i], end[j]);

                    double sum = 0.0;
                    for (int k = lk; k <= hk; k++)
                    {
                        sum += B[k][i - index[k]] * B[k][j - index[k]];
                    }
                    NTN[i][j] = sum;
                }
            }

            var Rk = new Vector3[m - 1];
            for (int k = 1; k <= m - 1; k++)
            {
                double n0 = Polynomials.OneBasisFunction(0, degree, knotVector, uk[k]);
                double np = Polynomials.OneBasisFunction(n, degree, knotVector, uk[k]);
                Rk[k - 1] =
                    throughPoints[k] - (float)n0 * throughPoints[0] - (float)np * throughPoints[m];
            }

            var right = new double[dim][];
            for (int i = 0; i < dim; i++)
            {
                Vector3 rSum = Vector3.Zero;
                int lk = start[i];
                int hk = end[i];

                for (int k = lk; k <= hk; k++)
                {
                    rSum += Rk[k] * (float)B[k][i - index[k]];
                }
                right[i] = new double[] { rSum.X, rSum.Y, rSum.Z };
            }

            try
            {
                var result = MathUtils.SolveLinearSystem(NTN, right);

                var controlPoints = new ControlPoint<Vector3>[controlPointsCount];
                controlPoints[0] = new ControlPoint<Vector3>(throughPoints[0], 1.0);
                controlPoints[n] = new ControlPoint<Vector3>(throughPoints[m], 1.0);

                for (int i = 1; i < n; i++)
                {
                    var row = result[i - 1];
                    controlPoints[i] = new ControlPoint<Vector3>(
                        new Vector3((float)row[0], (float)row[1], (float)row[2]),
                        1.0
                    );
                }

                curve = new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
                return true;
            }
            catch
            {
                return false;
            }
        }

        public static bool WeightedAndContrainedLeastSquaresApproximation(
            int degree,
            IReadOnlyList<Vector3> throughPoints,
            IReadOnlyList<double> throughPointWeights,
            IReadOnlyList<Vector3> tangents,
            IReadOnlyList<int> tangentIndices,
            IReadOnlyList<double> tangentWeights,
            int controlPointsCount,
            out NurbsCurve<Vector3> curve
        )
        {
            curve = default;
            Validate.Argument(degree > 0, nameof(degree), "Degree must be greater than zero.");
            int size = throughPoints.Count;
            Validate.Argument(
                size > degree,
                nameof(throughPoints),
                "ThroughPoints size must be greater than degree."
            );
            Validate.Argument(
                throughPointWeights.Count == size,
                nameof(throughPointWeights),
                "Weights size must be equal to throughPoints size."
            );

            int ru = -1;
            int rc = -1;
            int r = size - 1;
            for (int i = 0; i <= r; i++)
            {
                if (MathUtils.IsGreaterThan(throughPointWeights[i], 0.0))
                    ru++;
                else
                    rc++;
            }

            int su = -1;
            int sc = -1;
            int s = tangents.Count - 1;
            for (int j = 0; j <= s; j++)
            {
                if (MathUtils.IsGreaterThan(tangentWeights[j], 0.0))
                    su++;
                else
                    sc++;
            }

            int mu = ru + su + 1;
            int mc = rc + sc + 1;
            int n = controlPointsCount - 1;

            if (mc >= n || mc + n >= mu + 1)
                return false;

            var uk = Interpolation3D.GetChordParameterization(throughPoints);
            var knotVector = Interpolation.ComputeKnotVector(degree, controlPointsCount, uk);
            var controlPoints = new ControlPoint<Vector3>[controlPointsCount];

            int tangentIdx = 0;
            int mu2 = 0;
            int mc2 = 0;

            var N = new double[mu + 1][];
            for (int k = 0; k <= mu; k++)
                N[k] = new double[n + 1];

            var M = new double[mc + 1][];
            for (int k = 0; k <= mc; k++)
                M[k] = new double[n + 1];

            var S = new double[mu + 1][];
            var T = new double[mc + 1][];

            var W = new double[mu + 1];

            for (int i = 0; i <= r; i++)
            {
                int spanIndex = Polynomials.GetKnotSpanIndex(degree, knotVector, uk[i]);

                bool dflag = false;
                if (tangentIdx <= s)
                {
                    if (i == tangentIndices[tangentIdx])
                        dflag = true;
                }

                double[] funs0;
                double[] funs1 = null;

                if (!dflag)
                {
                    funs0 = Polynomials.BasisFunctions(spanIndex, degree, knotVector, uk[i]);
                }
                else
                {
                    var ders = Polynomials.BasisFunctionsDerivatives(
                        spanIndex,
                        degree,
                        1,
                        knotVector,
                        uk[i]
                    );
                    funs0 = ders[0];
                    funs1 = ders[1];
                }

                if (MathUtils.IsGreaterThan(throughPointWeights[i], 0.0))
                {
                    W[mu2] = throughPointWeights[i];
                    for (int z = 0; z < funs0.Length; z++)
                    {
                        N[mu2][spanIndex - degree + z] = funs0[z];
                    }
                    var sp = (float)throughPointWeights[i] * throughPoints[i];
                    S[mu2] = new double[] { sp.X, sp.Y, sp.Z };
                    mu2++;
                }
                else
                {
                    for (int z = 0; z < funs0.Length; z++)
                    {
                        M[mc2][spanIndex - degree + z] = funs0[z];
                    }
                    T[mc2] = new double[]
                    {
                        throughPoints[i].X,
                        throughPoints[i].Y,
                        throughPoints[i].Z,
                    };
                    mc2++;
                }

                if (dflag)
                {
                    if (MathUtils.IsGreaterThan(tangentWeights[tangentIdx], 0.0))
                    {
                        W[mu2] = tangentWeights[tangentIdx];
                        for (int z = 0; z < funs1.Length; z++)
                        {
                            N[mu2][spanIndex - degree + z] = funs1[z];
                        }
                        var sp = (float)tangentWeights[tangentIdx] * tangents[tangentIdx];
                        S[mu2] = new double[] { sp.X, sp.Y, sp.Z };
                        mu2++;
                    }
                    else
                    {
                        for (int z = 0; z < funs1.Length; z++)
                        {
                            M[mc2][spanIndex - degree + z] = funs1[z];
                        }
                        T[mc2] = new double[]
                        {
                            tangents[tangentIdx].X,
                            tangents[tangentIdx].Y,
                            tangents[tangentIdx].Z,
                        };
                        mc2++;
                    }
                    tangentIdx++;
                }
            }

            var tN = MathUtils.Transpose(N);
            var W_ = MathUtils.MakeDiagonal(W);
            var tNW = MathUtils.MatrixMultiply(tN, W_);
            var tNWN = MathUtils.MatrixMultiply(tNW, N);
            var tNWS = MathUtils.MatrixMultiply(tNW, S);

            double[][] resultMatrix;

            if (mc < 0)
            {
                resultMatrix = MathUtils.SolveLinearSystem(tNWN, tNWS);
            }
            else
            {
                if (!MathUtils.MakeInverse(tNWN, out var inv_tNWN))
                    return false;

                var tM = MathUtils.Transpose(M);

                var Minv_tNWN = MathUtils.MatrixMultiply(M, inv_tNWN);
                var Minv_tNWN_tM = MathUtils.MatrixMultiply(Minv_tNWN, tM);
                var Minv_tNWN_tNWS = MathUtils.MatrixMultiply(Minv_tNWN, tNWS);

                var RhsA = new double[T.Length][];
                for (int i = 0; i < T.Length; i++)
                {
                    RhsA[i] = new double[3];
                    RhsA[i][0] = Minv_tNWN_tNWS[i][0] - T[i][0];
                    RhsA[i][1] = Minv_tNWN_tNWS[i][1] - T[i][1];
                    RhsA[i][2] = Minv_tNWN_tNWS[i][2] - T[i][2];
                }

                var A = MathUtils.SolveLinearSystem(Minv_tNWN_tM, RhsA);
                var tMA = MathUtils.MatrixMultiply(tM, A);

                var RhsResult = new double[tNWS.Length][];
                for (int i = 0; i < tNWS.Length; i++)
                {
                    RhsResult[i] = new double[3];
                    RhsResult[i][0] = tNWS[i][0] - tMA[i][0];
                    RhsResult[i][1] = tNWS[i][1] - tMA[i][1];
                    RhsResult[i][2] = tNWS[i][2] - tMA[i][2];
                }

                resultMatrix = MathUtils.MatrixMultiply(inv_tNWN, RhsResult);
            }

            for (int i = 0; i < resultMatrix.Length; i++)
            {
                var pt = new Vector3(
                    (float)resultMatrix[i][0],
                    (float)resultMatrix[i][1],
                    (float)resultMatrix[i][2]
                );
                controlPoints[i] = new ControlPoint<Vector3>(pt, 1.0);
            }

            curve = new NurbsCurve<Vector3>(degree, controlPoints, knotVector);
            return true;
        }

        private static double ComputeRemoveKnotErrorBound(
            NurbsCurve<Vector3> curve,
            int removalIndex
        )
        {
            int degree = curve.Degree;
            var knotVector = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Range(removalIndex, 0, knotVector.Count - 1, nameof(removalIndex));

            int ord = degree + 1;
            int r = removalIndex;
            double u = knotVector[r];
            int s = Polynomials.GetKnotMultiplicity(knotVector, u);
            int last = r - s;
            int first = r - degree;
            int off = first - 1;

            var temp = new Vector4[knotVector.Count];

            var cpOff = controlPoints[off];
            // Homogeneous: x*w, y*w, z*w, w
            temp[0] = new Vector4(cpOff.Value * (float)cpOff.Weight, (float)cpOff.Weight);

            var cpLast = controlPoints[last + 1];
            temp[last + 1 - off] = new Vector4(
                cpLast.Value * (float)cpLast.Weight,
                (float)cpLast.Weight
            );

            int i = first;
            int j = last;
            int ii = 1;
            int jj = last - off;

            while (j - i > 0)
            {
                double alfi = (u - knotVector[i]) / (knotVector[i + ord] - knotVector[i]);
                double alfj = (u - knotVector[j]) / (knotVector[j + ord] - knotVector[j]);

                var pI = controlPoints[i];
                var vI = new Vector4(pI.Value * (float)pI.Weight, (float)pI.Weight);
                temp[ii] = (vI - (float)(1.0 - alfi) * temp[ii - 1]) / (float)alfi;

                var pJ = controlPoints[j];
                var vJ = new Vector4(pJ.Value * (float)pJ.Weight, (float)pJ.Weight);
                temp[jj] = (vJ - (float)alfj * temp[jj + 1]) / (float)(1.0 - alfj);

                i++;
                ii++;
                j--;
                jj--;
            }

            if (j - i < 0)
            {
                return (temp[ii - 1] - temp[jj + 1]).Length();
            }
            else
            {
                double alfi = (u - knotVector[i]) / (knotVector[i + ord] - knotVector[i]);
                var cpI = controlPoints[i];
                var vI = new Vector4(cpI.Value * (float)cpI.Weight, (float)cpI.Weight);
                return (
                    vI - ((float)alfi * temp[ii + 1] + (float)(1.0 - alfi) * temp[ii - 1])
                ).Length();
            }
        }

        public static bool RemoveKnotsByGivenBound(
            NurbsCurve<Vector3> curve,
            IReadOnlyList<double> parameters,
            List<double> errors,
            double maxError,
            out NurbsCurve<Vector3> result
        )
        {
            result = default;
            int degree = curve.Degree;
            // Create mutable copies
            var tempU = new List<double>(curve.Knots);
            var tempCP = new List<ControlPoint<Vector3>>(curve.ControlPoints);
            var currentCurve = new NurbsCurve<Vector3>(degree, tempCP, tempU);

            Validate.Argument(
                parameters.Count > 0,
                nameof(parameters),
                "Params size must be greater than zero."
            );
            Validate.Argument(
                parameters.Count == errors.Count,
                nameof(errors),
                "Errors size must be equal to params size."
            );
            Validate.Argument(
                MathUtils.IsGreaterThan(maxError, 0.0),
                nameof(maxError),
                "Maxerror must be greater than zero."
            );

            int knotSize = tempU.Count;
            var Br = new List<double>(new double[knotSize]);
            var S = new List<int>(new int[knotSize]);
            var Nl = new List<int>(new int[knotSize]);
            var Nr = new List<int>(new int[knotSize]);

            for (int k = 0; k < knotSize; k++)
            {
                Br[k] = Constants.MaxDistance;
                Nr[k] = parameters.Count - 1;
            }

            var uk = parameters;
            int ukSize = uk.Count;
            var NewError = new double[ukSize];
            var temp = new double[ukSize];

            int s = 1;
            int controlPointsSize = tempCP.Count;
            int n = controlPointsSize - 1;

            for (int i = degree + 1; i < controlPointsSize; i++)
            {
                if (MathUtils.IsLessThan(tempU[i], tempU[i + 1]))
                {
                    Br[i] = ComputeRemoveKnotErrorBound(currentCurve, i);
                    S[i] = Polynomials.GetKnotMultiplicity(tempU, tempU[i]);
                    s = 1;
                }
                else
                {
                    Br[i] = Constants.MaxDistance;
                    S[i] = 1;
                    s++;
                }
            }

            Nl[0] = 0;
            for (int i = 0; i < ukSize; i++)
            {
                int spanIndex = Polynomials.GetKnotSpanIndex(degree, tempU, uk[i]);
                if (Nl[spanIndex] == 0 && spanIndex != 0)
                {
                    Nl[spanIndex] = i;
                }
                if (i + 1 < ukSize)
                {
                    Nr[spanIndex] = i + 1;
                }
            }

            while (true)
            {
                double minStandard = Constants.MaxDistance;
                int BrMinIndex = 0;

                for (int i = 0; i < Br.Count; i++)
                {
                    if (MathUtils.IsLessThan(Br[i], minStandard))
                    {
                        BrMinIndex = i;
                        minStandard = Br[i];
                    }
                }

                double BrMin = Br[BrMinIndex];

                if (MathUtils.IsAlmostEqualTo(BrMin, Constants.MaxDistance))
                {
                    break;
                }

                int r = BrMinIndex;
                s = S[BrMinIndex];

                int Rstart = Math.Max(r - degree, degree + 1);
                int Rend = Math.Min(r + degree - S[r + degree] + 1, n);
                Rstart = Nl[Rstart];
                Rend = Nr[Rend];

                bool removable = true;
                for (int i = Rstart; i <= Rend; i++)
                {
                    double a;
                    if ((degree + s) % 2 != 0)
                    {
                        double u = uk[i];
                        int k = (degree + s + 1) / 2;
                        a = tempU[r] - tempU[r - k + 1];
                        a /= tempU[r - k + degree + 2] - tempU[r - k + 1];
                        NewError[i] =
                            (1.0 - a)
                            * Br[r]
                            * Polynomials.OneBasisFunction(r - k + 1, degree, tempU, u);
                    }
                    else
                    {
                        double u = uk[i];
                        int k = (degree + s) / 2;
                        NewError[i] = Br[r] * Polynomials.OneBasisFunction(r - k, degree, tempU, u);
                    }
                    temp[i] = NewError[i] + errors[i];
                    if (MathUtils.IsGreaterThan(temp[i], maxError))
                    {
                        removable = false;
                        Br[r] = Constants.MaxDistance;
                        break;
                    }
                }

                if (removable)
                {
                    if (RemoveKnot(currentCurve, tempU[r], 1, out var newtc))
                    {
                        currentCurve = newtc;
                        tempU = new List<double>(newtc.Knots);
                        tempCP = new List<ControlPoint<Vector3>>(newtc.ControlPoints);

                        controlPointsSize = tempCP.Count;
                        n = controlPointsSize - 1;

                        for (int i = Rstart; i <= Rend; i++)
                        {
                            errors[i] = temp[i];
                        }

                        if (controlPointsSize <= degree + 1)
                        {
                            break;
                        }

                        Rstart = Nl[r - degree - 1];
                        Rend = Nr[r - S[r]];

                        int spanIndex = 0;
                        int oldspanIndex = -1;
                        for (int k = Rstart; k <= Rend; k++)
                        {
                            spanIndex = Polynomials.GetKnotSpanIndex(degree, tempU, uk[k]);
                            if (spanIndex != oldspanIndex)
                            {
                                Nl[spanIndex] = k;
                            }
                            if (k + 1 < ukSize)
                            {
                                Nr[spanIndex] = k + 1;
                            }
                            oldspanIndex = spanIndex;
                        }

                        for (int k = r - S[r] + 1; k < Nl.Count - 1; k++)
                        {
                            Nl[k] = Nl[k + 1];
                            Nr[k] = Nr[k + 1];
                        }
                        Nl.RemoveAt(Nl.Count - 1);
                        Nr.RemoveAt(Nr.Count - 1);

                        Rstart = Math.Max(r - degree, degree + 1);
                        Rend = Math.Min(r + degree - S[r] + 1, controlPointsSize);
                        s = S[Rstart];

                        for (int i = Rstart; i <= Rend; i++)
                        {
                            if (MathUtils.IsLessThan(tempU[i], tempU[i + 1]))
                            {
                                Br[i] = ComputeRemoveKnotErrorBound(currentCurve, i);
                                S[i] = s;
                                s = 1;
                            }
                            else
                            {
                                Br[i] = Constants.MaxDistance;
                                S[i] = 1;
                                s++;
                            }
                        }

                        for (int i = Rend + 1; i < Br.Count - 1; i++)
                        {
                            Br[i] = Br[i + 1];
                            S[i] = S[i + 1];
                        }
                        Br.RemoveAt(Br.Count - 1);
                        S.RemoveAt(S.Count - 1);
                    }
                    else
                    {
                        Br[r] = Constants.MaxDistance;
                    }
                }
                else
                {
                    Br[r] = Constants.MaxDistance;
                }
            }

            result = currentCurve;
            return true;
        }

        public static bool GlobalApproximationByErrorBound(
            int degree,
            IReadOnlyList<Vector3> throughPoints,
            double maxError,
            out NurbsCurve<Vector3> curve
        )
        {
            curve = default;
            Validate.Argument(degree > 0, nameof(degree), "Degree must be greater than zero.");
            int size = throughPoints.Count;
            Validate.Argument(
                size > degree,
                nameof(throughPoints),
                "ThroughPoints size must be greater than degree."
            );
            Validate.Argument(
                MathUtils.IsGreaterThan(maxError, 0.0),
                nameof(maxError),
                "Maxerror must be greater than zero."
            );

            var uk = Interpolation.GetChordParameterization(throughPoints);
            var errors = new List<double>(new double[size]);

            var controlPoints = new ControlPoint<Vector3>[size];
            for (int i = 0; i < size; i++)
            {
                controlPoints[i] = new ControlPoint<Vector3>(throughPoints[i], 1.0);
            }

            var knotVector = new double[size + 2];
            for (int i = 0; i < size; i++)
            {
                knotVector[i + 1] = uk[i];
            }
            knotVector[0] = 0.0;
            knotVector[size + 1] = 1.0;

            var tc = new NurbsCurve<Vector3>(1, controlPoints, knotVector);
            var newtc = ElevateDegree(tc, degree - 1);

            return RemoveKnotsByGivenBound(newtc, uk, errors, maxError, out curve);
        }

        public static bool FitWithConic(
            IReadOnlyList<Vector3> throughPoints,
            int startPointIndex,
            int endPointIndex,
            Vector3 startTangent,
            Vector3 endTangent,
            double maxError,
            List<ControlPoint<Vector3>> middleControlPoints
        )
        {
            Validate.Argument(
                throughPoints.Count > 0,
                nameof(throughPoints),
                "ThroughPoints size must be greater than zero."
            );
            Validate.Range(startPointIndex, 0, throughPoints.Count - 1, nameof(startPointIndex));
            Validate.Range(
                endPointIndex,
                startPointIndex + 1,
                throughPoints.Count - 1,
                nameof(endPointIndex)
            );
            Validate.Argument(
                !MathUtils.IsZero(startTangent),
                nameof(startTangent),
                "StartTangent must not be zero vector."
            );
            Validate.Argument(
                !MathUtils.IsZero(endTangent),
                nameof(endTangent),
                "EndTangent must not be zero vector."
            );

            var startPoint = throughPoints[startPointIndex];
            var endPoint = throughPoints[endPointIndex];

            if (endPointIndex - startPointIndex == 1)
            {
                var computed = BezierCurve<Vector3>.ComputeMiddleControlPointsOnQuadraticCurve(
                    startPoint,
                    startTangent,
                    endPoint,
                    endTangent
                );
                if (computed != null)
                {
                    foreach (var cp in computed)
                    {
                        middleControlPoints.Add(cp);
                    }
                    return true;
                }
                return false;
            }

            var type = Intersection.ComputeRays(
                startPoint,
                startTangent,
                endPoint,
                endTangent,
                out double alf1,
                out double alf2,
                out Vector3 R
            );

            if (type == CurveCurveIntersectionType.Coincident)
            {
                middleControlPoints.Add(
                    new ControlPoint<Vector3>((startPoint + endPoint) * 0.5f, 1.0)
                );
                return true;
            }
            else if (
                type == CurveCurveIntersectionType.Skew
                || type == CurveCurveIntersectionType.Parallel
            )
            {
                return false;
            }

            if (MathUtils.IsLessThanOrEqual(alf1, 0.0) || MathUtils.IsGreaterThanOrEqual(alf2, 0.0))
            {
                return false;
            }

            double s = 0.0;
            Vector3 V = endPoint - startPoint;

            for (int i = startPointIndex + 1; i <= endPointIndex - 1; i++)
            {
                Vector3 V1 = throughPoints[i] - R;
                type = Intersection.ComputeRays(
                    startPoint,
                    V,
                    R,
                    V1,
                    out double a1,
                    out double a2,
                    out Vector3 dummy
                );

                if (
                    type != CurveCurveIntersectionType.Intersecting
                    || MathUtils.IsLessThanOrEqual(a1, 0.0)
                    || MathUtils.IsGreaterThanOrEqual(a1, 1.0)
                    || MathUtils.IsLessThanOrEqual(a2, 0.0)
                )
                {
                    return false;
                }

                if (
                    CreateOneConicArc(
                        startPoint,
                        V,
                        R,
                        V1,
                        throughPoints[i],
                        out dummy,
                        out double wi
                    )
                )
                {
                    s = s + wi / (1.0 + wi);
                }
            }

            s = s / (endPointIndex - startPointIndex - 1);
            double w = s / (1.0 - s);

            var controlPoints = new ControlPoint<Vector3>[]
            {
                new ControlPoint<Vector3>(startPoint, 1.0),
                new ControlPoint<Vector3>(R, w),
                new ControlPoint<Vector3>(endPoint, 1.0),
            };
            var knotVectors = new double[] { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };

            var tc = new NurbsCurve<Vector3>(2, controlPoints, knotVectors);

            for (int k = startPointIndex + 1; k <= endPointIndex - 1; k++)
            {
                var tp = throughPoints[k];
                double param = GetParamOnCurve(tc, tp);
                Vector3 point = GetPointOnCurve(tc, param);
                if (MathUtils.IsGreaterThan((tp - point).Length(), maxError))
                {
                    return false;
                }
            }

            middleControlPoints.Add(new ControlPoint<Vector3>(R, w));
            return true;
        }

        public static bool FitWithCubic(
            IReadOnlyList<Vector3> throughPoints,
            int startPointIndex,
            int endPointIndex,
            Vector3 startTangent,
            Vector3 endTangent,
            double maxError,
            List<ControlPoint<Vector3>> middleControlPoints
        )
        {
            Validate.Argument(
                throughPoints.Count >= 3,
                nameof(throughPoints),
                "ThroughPoints size must be greater than 2."
            );
            Validate.Range(startPointIndex, 0, throughPoints.Count - 1, nameof(startPointIndex));
            Validate.Range(
                endPointIndex,
                startPointIndex + 1,
                throughPoints.Count - 1,
                nameof(endPointIndex)
            );
            Validate.Argument(
                !MathUtils.IsZero(startTangent),
                nameof(startTangent),
                "StartTangent must not be zero vector."
            );
            Validate.Argument(
                !MathUtils.IsZero(endTangent),
                nameof(endTangent),
                "EndTangent must not be zero vector."
            );

            var startPoint = throughPoints[startPointIndex];
            var endPoint = throughPoints[endPointIndex];
            int size = throughPoints.Count;

            if (endPointIndex - startPointIndex == 1)
            {
                if (!Interpolation.TryComputeTangents(throughPoints, out var tangents))
                    return false;

                Vector3 dks;
                if (startPointIndex == 0)
                {
                    dks = tangents[startPointIndex];
                }
                else
                {
                    double d1 = Vector3.Distance(endPoint, startPoint);
                    double d2 = Vector3.Distance(startPoint, throughPoints[startPointIndex - 1]);
                    dks = (float)(d1 / d2) * startTangent;
                }

                Vector3 dke;
                if (endPointIndex == size - 1)
                {
                    dke = tangents[endPointIndex];
                }
                else
                {
                    double d1 = Vector3.Distance(throughPoints[endPointIndex + 1], endPoint);
                    double d2 = Vector3.Distance(endPoint, startPoint);
                    dke = (float)(d1 / d2) * endTangent;
                }

                double alpha = dks.Length() / 3.0;
                double beta = -dke.Length() / 3.0;

                var p1 = startPoint + (float)alpha * startTangent;
                var p2 = endPoint + (float)beta * endTangent;

                middleControlPoints.Add(new ControlPoint<Vector3>(p1, 1.0));
                middleControlPoints.Add(new ControlPoint<Vector3>(p2, 1.0));
                return true;
            }

            int dk = endPointIndex - startPointIndex;
            bool isLine = true;
            Vector3 tempStandard = Vector3.Zero;

            for (int i = startPointIndex; i <= endPointIndex; i++)
            {
                var tp = throughPoints[i];
                if (MathUtils.IsAlmostEqualTo(tp, startPoint))
                    continue;

                var direction = tp - startPoint;
                direction.Normalize();

                if (MathUtils.IsZero(tempStandard))
                {
                    tempStandard = direction;
                }
                else
                {
                    if (!MathUtils.IsAlmostEqualTo(tempStandard, direction))
                    {
                        isLine = false;
                        break;
                    }
                }
            }

            if (isLine)
            {
                var p1 = (2.0f * startPoint + endPoint) / 3.0f;
                var p2 = (startPoint + 2.0f * endPoint) / 3.0f;
                middleControlPoints.Add(new ControlPoint<Vector3>(p1, 1.0));
                middleControlPoints.Add(new ControlPoint<Vector3>(p2, 1.0));
                return true;
            }

            var newThroughPoints = new List<Vector3>(dk + 1);
            for (int i = startPointIndex; i <= endPointIndex; i++)
                newThroughPoints.Add(throughPoints[i]);

            var uh = Interpolation.GetChordParameterization(newThroughPoints);
            var alphak = new double[dk + 1];
            var betak = new double[dk + 1];
            bool possible = true;

            for (int k = 1; k < dk; k++)
            {
                double u = uh[k];
                double s = 1.0 - u;

                double b0 = s * s * s;
                double b1 = 3 * s * s * u;
                double b2 = 3 * s * u * u;
                double b3 = u * u * u;

                var termP0 = (float)(b0 + b1) * startPoint;
                var termP3 = (float)(b2 + b3) * endPoint;
                var rhs = throughPoints[startPointIndex + k] - termP0 - termP3;

                var vecA = (float)b1 * startTangent;
                var vecB = (float)b2 * endTangent;

                // Least squares for 3x2: [A B] [alpha; beta] = RHS
                double aa = Vector3.Dot(vecA, vecA);
                double ab = Vector3.Dot(vecA, vecB);
                double bb = Vector3.Dot(vecB, vecB);
                double ar = Vector3.Dot(vecA, rhs);
                double br = Vector3.Dot(vecB, rhs);

                double det = aa * bb - ab * ab;
                if (Math.Abs(det) < MathUtils.Epsilon)
                {
                    possible = false;
                    break;
                }

                double ak = (ar * bb - br * ab) / det;
                double bk = (aa * br - ab * ar) / det;

                if (ak > 0.0 && bk < 0.0)
                {
                    alphak[k] = ak;
                    betak[k] = bk;
                }
                else
                {
                    possible = false;
                    break;
                }
            }

            if (!possible)
                return false;

            double alphaAvg = 0.0;
            double betaAvg = 0.0;
            for (int k = 1; k < dk; k++)
            {
                alphaAvg += alphak[k];
                betaAvg += betak[k];
            }
            alphaAvg /= (dk - 1);
            betaAvg /= (dk - 1);

            var finalP1 = startPoint + (float)alphaAvg * startTangent;
            var finalP2 = endPoint + (float)betaAvg * endTangent;

            var controlPoints = new ControlPoint<Vector3>[]
            {
                new ControlPoint<Vector3>(startPoint, 1.0),
                new ControlPoint<Vector3>(finalP1, 1.0),
                new ControlPoint<Vector3>(finalP2, 1.0),
                new ControlPoint<Vector3>(endPoint, 1.0),
            };
            var curve = new BezierCurve<Vector3>(controlPoints);

            for (int k = 1; k < dk; k++)
            {
                double u = uh[k];
                var p = BezierCurve3DHelper.GetPointOnCurveByBernstein(curve, u);
                if (
                    MathUtils.IsGreaterThan(
                        Vector3.Distance(throughPoints[startPointIndex + k], p),
                        maxError
                    )
                )
                {
                    return false;
                }
            }

            middleControlPoints.Add(new ControlPoint<Vector3>(finalP1, 1.0));
            middleControlPoints.Add(new ControlPoint<Vector3>(finalP2, 1.0));
            return true;
        }

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

        private static void TessellateCore(
            NurbsCurve<Vector3> curve,
            double start,
            double end,
            List<double> parameters,
            int depth
        )
        {
            const int MaxDepth = 20;

            if (depth >= MaxDepth)
            {
                parameters.Add((start + end) * 0.5);
                return;
            }

            var p0 = GetPointOnCurve(curve, start);
            var p2 = GetPointOnCurve(curve, end);
            double mid = (start + end) * 0.5;
            var p1 = GetPointOnCurve(curve, mid);

            var chordMid = (p0 + p2) * 0.5f;
            double deviation = Vector3.Distance(p1, chordMid);

            if (MathUtils.IsAlmostEqualTo(deviation, 0.0))
            {
                parameters.Add(mid);
            }
            else
            {
                TessellateCore(curve, start, mid, parameters, depth + 1);
                TessellateCore(curve, mid, end, parameters, depth + 1);
            }
        }

        public static IReadOnlyList<Vector3> Tessellate(NurbsCurve<Vector3> curve)
        {
            if (curve.Degree == 1)
            {
                var result = new Vector3[curve.ControlPoints.Count];
                for (int i = 0; i < curve.ControlPoints.Count; i++)
                {
                    result[i] = curve.ControlPoints[i].Value;
                }
                return result;
            }

            var uniqueKnots = new List<double>();
            var knots = curve.Knots;
            if (knots.Count > 0)
            {
                uniqueKnots.Add(knots[0]);
                for (int i = 1; i < knots.Count; i++)
                {
                    if (!MathUtils.IsAlmostEqualTo(knots[i], uniqueKnots[uniqueKnots.Count - 1]))
                    {
                        uniqueKnots.Add(knots[i]);
                    }
                }
            }

            if (uniqueKnots.Count < 2)
            {
                var pt = GetPointOnCurve(curve, uniqueKnots[0]);
                return new Vector3[] { pt };
            }

            var parameters = new List<double>();
            double u_start = uniqueKnots[0];
            double u_end = uniqueKnots[uniqueKnots.Count - 1];

            parameters.Add(u_start);

            for (int i = 0; i < uniqueKnots.Count - 1; ++i)
            {
                double u0 = uniqueKnots[i];
                double u1 = uniqueKnots[i + 1];

                if (MathUtils.IsAlmostEqualTo(u1, u0))
                {
                    continue;
                }

                var internalParams = new List<double>();
                TessellateCore(curve, u0, u1, internalParams, 0);

                foreach (double t in internalParams)
                {
                    if (parameters.Count == 0)
                        continue;
                    if (
                        t > u0
                        && t < u1
                        && (MathUtils.IsGreaterThan(t, parameters[parameters.Count - 1]))
                    )
                    {
                        parameters.Add(t);
                    }
                }

                if (MathUtils.IsGreaterThan(u1, parameters[parameters.Count - 1]))
                {
                    parameters.Add(u1);
                }
            }

            if (
                parameters.Count == 0
                || !MathUtils.IsAlmostEqualTo(parameters[parameters.Count - 1], u_end)
            )
            {
                parameters.Add(u_end);
            }

            var points = new Vector3[parameters.Count];
            for (int i = 0; i < parameters.Count; i++)
            {
                points[i] = GetPointOnCurve(curve, parameters[i]);
            }
            return points;
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

        public static bool ControlPointReposition(
            NurbsCurve<Vector3> curve,
            double parameter,
            int moveIndex,
            Vector3 moveDirection,
            double moveDistance,
            out NurbsCurve<Vector3> result
        )
        {
            result = default;
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Range(parameter, knots[0], knots[knots.Count - 1], nameof(parameter));
            Validate.Range(moveIndex, 0, controlPoints.Count - 1, nameof(moveIndex));
            Validate.Argument(
                !MathUtils.IsZero(moveDirection),
                nameof(moveDirection),
                "MoveDirection must not be zero vector."
            );
            Validate.Argument(
                !MathUtils.IsAlmostEqualTo(moveDistance, 0.0),
                nameof(moveDistance),
                "MoveDistance must not be zero."
            );

            int spanIndex = Polynomials.GetKnotSpanIndex(degree, knots, parameter);

            if (moveIndex < spanIndex - degree || moveIndex > spanIndex)
            {
                return false;
            }

            var basis = Polynomials.BasisFunctions(spanIndex, degree, knots, parameter);

            double den = 0.0;
            for (int i = 0; i <= degree; i++)
            {
                den += basis[i] * controlPoints[spanIndex - degree + i].Weight;
            }

            if (MathUtils.IsZero(den))
                return false;

            double num = basis[moveIndex - (spanIndex - degree)] * controlPoints[moveIndex].Weight;
            double Rkp = num / den;

            if (MathUtils.IsLessThan(Rkp, 0.0) || MathUtils.IsZero(Rkp))
            {
                return false;
            }

            var updatedControlPoints = new List<ControlPoint<Vector3>>(controlPoints);
            Vector3 movePoint = controlPoints[moveIndex].Value;
            double alpha = moveDistance / (moveDirection.Length() * Rkp);
            Vector3 newPoint = movePoint + (float)alpha * moveDirection;

            updatedControlPoints[moveIndex] = new ControlPoint<Vector3>(
                newPoint,
                controlPoints[moveIndex].Weight
            );

            result = new NurbsCurve<Vector3>(degree, updatedControlPoints, knots);
            return true;
        }

        public static bool NeighborWeightsModification(
            NurbsCurve<Vector3> curve,
            double parameter,
            int moveIndex,
            double moveDistance,
            double scale,
            out NurbsCurve<Vector3> result
        )
        {
            result = default;
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Range(parameter, knots[0], knots[knots.Count - 1], nameof(parameter));
            Validate.Range(moveIndex, 0, controlPoints.Count - 2, nameof(moveIndex));
            Validate.Argument(
                !MathUtils.IsAlmostEqualTo(moveDistance, 0.0),
                nameof(moveDistance),
                "MoveDistance must not be zero."
            );
            Validate.Argument(
                !MathUtils.IsAlmostEqualTo(scale, 0.0),
                nameof(scale),
                "Scale must not be zero."
            );

            var tempControlPoints = new List<ControlPoint<Vector3>>(controlPoints);
            var movePoint1 = tempControlPoints[moveIndex].Value;
            var movePoint2 = tempControlPoints[moveIndex + 1].Value;

            tempControlPoints[moveIndex] = new ControlPoint<Vector3>(movePoint1, 0.0);
            tempControlPoints[moveIndex + 1] = new ControlPoint<Vector3>(movePoint2, 0.0);

            var tc = new NurbsCurve<Vector3>(degree, tempControlPoints, knots);
            var R = GetPointOnCurve(tc, parameter);

            var controlLeg = movePoint1 - movePoint2;
            var controlLegLength = Vector3.Distance(movePoint1, movePoint2);

            var P = GetPointOnCurve(curve, parameter);
            var direction = R - P;

            var type = Intersection.ComputeRays(
                movePoint1,
                controlLeg,
                R,
                direction,
                out _,
                out _,
                out Vector3 Q
            );

            if (type != CurveCurveIntersectionType.Intersecting)
                return false;

            var pkq = Q - movePoint1;
            var pk1q = Q - movePoint2;

            double RQ = Vector3.Distance(Q, R);
            double RP = Vector3.Distance(P, R);

            if (MathUtils.IsZero(RP) || MathUtils.IsZero(RQ) || MathUtils.IsZero(controlLegLength))
                return false;

            double Rtarget = RP + moveDistance;
            double qRP = RP / RQ;
            double qRtarget = Rtarget / RQ;

            var A = movePoint1 + (float)qRP * pkq;
            var B = movePoint2 + (float)qRP * pk1q;
            var C = movePoint1 + (float)qRtarget * pkq;
            var D = movePoint2 + (float)qRtarget * pk1q;

            double ak = Vector3.Distance(B, movePoint2) / controlLegLength;
            double ak1 = Vector3.Distance(A, movePoint1) / controlLegLength;
            double abk = Vector3.Distance(D, movePoint2) / controlLegLength;
            double abk1 = Vector3.Distance(C, movePoint1) / controlLegLength;

            if (
                MathUtils.IsLessThan(Math.Abs(ak), 0.0)
                || MathUtils.IsLessThan(Math.Abs(ak1), 0.0)
                || MathUtils.IsLessThan(Math.Abs(abk), 0.0)
                || MathUtils.IsLessThan(Math.Abs(abk1), 0.0)
            )
            {
                return false;
            }

            double alpha = 1.0 - ak - ak1;
            double beta = 1.0 - abk - abk1;

            if (
                MathUtils.IsZero(ak)
                || MathUtils.IsZero(abk)
                || MathUtils.IsZero(ak1)
                || MathUtils.IsZero(abk1)
            )
                return false;

            double betak = (alpha / ak) / (beta / abk);
            double betak1 = (alpha / ak1) / (beta / abk1);

            var updatedControlPoints = new List<ControlPoint<Vector3>>(controlPoints);
            updatedControlPoints[moveIndex] = new ControlPoint<Vector3>(
                movePoint1,
                controlPoints[moveIndex].Weight * betak
            );
            updatedControlPoints[moveIndex + 1] = new ControlPoint<Vector3>(
                movePoint2,
                controlPoints[moveIndex + 1].Weight * betak1
            );

            result = new NurbsCurve<Vector3>(degree, updatedControlPoints, knots);
            return true;
        }

        public static NurbsCurve<Vector3> Warping(
            NurbsCurve<Vector3> curve,
            IReadOnlyList<double> warpShape,
            double warpDistance,
            Vector3 planeNormal,
            double startParameter,
            double endParameter
        )
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Argument(
                controlPoints.Count == warpShape.Count,
                nameof(warpShape),
                "WarpShape size must be equal to control points size."
            );
            Validate.Argument(
                !MathUtils.IsAlmostEqualTo(warpDistance, 0.0),
                nameof(warpDistance),
                "WarpDistance must not be zero."
            );
            Validate.Argument(
                !MathUtils.IsZero(planeNormal),
                nameof(planeNormal),
                "PlaneNormal must not be zero vector."
            );
            Validate.Range(
                startParameter,
                knots[0],
                knots[knots.Count - 1],
                nameof(startParameter)
            );
            Validate.Range(
                endParameter,
                startParameter,
                knots[knots.Count - 1],
                nameof(endParameter)
            );
            Validate.Argument(
                MathUtils.IsGreaterThan(endParameter, startParameter),
                nameof(endParameter),
                "EndParameter must be greater than startParameter."
            );

            double halfParameter = 0.5 * (startParameter + endParameter);
            var derivatives = ComputeRationalCurveDerivatives(curve, 1, halfParameter);
            var tangent = derivatives[1];
            var normal = Vector3.Cross(tangent, planeNormal);

            var W = MathUtils.IsGreaterThan(warpDistance, 0.0) ? normal : -normal;
            W = Vector3.Normalize(W);

            var resultControlPoints = new ControlPoint<Vector3>[controlPoints.Count];
            double absWarpDistance = Math.Abs(warpDistance);

            for (int i = 0; i < controlPoints.Count; i++)
            {
                var cp = controlPoints[i];
                var currentPoint = cp.Value;
                var newPoint = currentPoint + (float)(warpShape[i] * absWarpDistance) * W;
                resultControlPoints[i] = new ControlPoint<Vector3>(newPoint, cp.Weight);
            }

            return new NurbsCurve<Vector3>(degree, resultControlPoints, knots);
        }

        public static bool Flattening(
            NurbsCurve<Vector3> curve,
            Vector3 lineStart,
            Vector3 lineEnd,
            double startParameter,
            double endParameter,
            out NurbsCurve<Vector3> result
        )
        {
            int degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Argument(
                !MathUtils.IsAlmostEqualTo(lineStart, lineEnd),
                nameof(lineEnd),
                "Line end point must not be equal to line start point."
            );
            Validate.Range(
                startParameter,
                knots[0],
                knots[knots.Count - 1],
                nameof(startParameter)
            );
            Validate.Range(
                endParameter,
                startParameter,
                knots[knots.Count - 1],
                nameof(endParameter)
            );
            Validate.Argument(
                endParameter > startParameter,
                nameof(endParameter),
                "EndParameter must be greater than StartParameter."
            );

            int spanMinIndex = Polynomials.GetKnotSpanIndex(degree, knots, startParameter);
            int spanMaxIndex = Polynomials.GetKnotSpanIndex(degree, knots, endParameter);

            var selectedControlPoints = new Dictionary<int, Vector3>();
            for (int i = spanMinIndex; i <= spanMaxIndex - degree - 1; i++)
            {
                if (i >= 0 && i < controlPoints.Count)
                {
                    selectedControlPoints.Add(i, controlPoints[i].Value);
                }
            }

            int projectCount = 0;
            var updatedControlPoints = new List<ControlPoint<Vector3>>(controlPoints);

            Vector3 lineVec = lineEnd - lineStart;
            double lineLenSq = lineVec.LengthSquared();

            foreach (var kvp in selectedControlPoints)
            {
                int index = kvp.Key;
                Vector3 current = kvp.Value;

                double t = Vector3.Dot(current - lineStart, lineVec) / (float)lineLenSq;
                if (t >= 0.0 && t <= 1.0)
                {
                    projectCount++;
                    var projected = lineStart + lineVec * (float)t;
                    updatedControlPoints[index] = new ControlPoint<Vector3>(
                        projected,
                        updatedControlPoints[index].Weight
                    );
                }
            }

            result = new NurbsCurve<Vector3>(degree, updatedControlPoints, knots);
            return projectCount >= degree + 1;
        }

        public static NurbsCurve<Vector3> Bending(
            NurbsCurve<Vector3> curve,
            double startParameter,
            double endParameter,
            Vector3 bendCenter,
            double radius,
            double crossRatio
        )
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Range(
                startParameter,
                knots[0],
                knots[knots.Count - 1],
                nameof(startParameter)
            );
            Validate.Range(
                endParameter,
                startParameter,
                knots[knots.Count - 1],
                nameof(endParameter)
            );

            int spanMinIndex = Polynomials.GetKnotSpanIndex(degree, knots, startParameter);
            int spanMaxIndex = Polynomials.GetKnotSpanIndex(degree, knots, endParameter);

            var updatedControlPoints = new List<ControlPoint<Vector3>>(controlPoints);

            for (int i = spanMinIndex; i <= spanMaxIndex - degree - 1; i++)
            {
                if (i < 0 || i >= updatedControlPoints.Count)
                    continue;

                var currentCP = updatedControlPoints[i];
                Vector3 current = currentCP.Value;

                var diff = current - bendCenter;

                if (diff.LengthSquared() < Constants.DoubleEpsilon)
                    continue;

                var pointOnBendCurve = bendCenter + Vector3.Normalize(diff) * (float)radius;

                double distCurrent = Vector3.Distance(bendCenter, current);
                double distBend = Vector3.Distance(bendCenter, pointOnBendCurve);

                double si = distBend / distCurrent;
                double ti = (crossRatio * si) / (1.0 + (crossRatio - 1.0) * si);

                Vector3 project = bendCenter + diff * (float)ti;

                updatedControlPoints[i] = new ControlPoint<Vector3>(project, currentCP.Weight);
            }

            return new NurbsCurve<Vector3>(degree, updatedControlPoints, knots);
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

        public static NurbsCurve<Vector3> ToClampCurve(NurbsCurve<Vector3> curve)
        {
            int degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            int m = knots.Count - 1;
            double up = knots[degree];
            double ump = knots[m - degree];

            int upMulti = Polynomials.GetKnotMultiplicity(knots, up);
            int umpMulti = Polynomials.GetKnotMultiplicity(knots, ump);

            int t1 = degree + 1 - upMulti;
            int t2 = degree + 1 - umpMulti;

            NurbsCurve<Vector3> tc = curve;
            if (t1 > 0)
            {
                InsertKnot(tc, up, t1, out tc);
            }
            if (t2 > 0)
            {
                InsertKnot(tc, ump, t2, out tc);
            }

            var kv = tc.Knots;
            var cps = tc.ControlPoints;
            int knotsCount = knots.Count;
            int cpCount = controlPoints.Count;

            var newKnots = new double[knotsCount];
            var newCPs = new ControlPoint<Vector3>[cpCount];

            for (int i = 0; i < knotsCount; i++)
            {
                newKnots[i] = kv[i + degree];
            }
            for (int i = 0; i < cpCount; i++)
            {
                newCPs[i] = cps[i + degree];
            }

            var result = new NurbsCurve<Vector3>(degree, newCPs, newKnots);

            bool isClosed = IsClosed(curve);
            bool isNowClosed = IsClosed(result);
            if (isClosed && !isNowClosed)
            {
                newCPs[newCPs.Length - 1] = newCPs[0];
                result = new NurbsCurve<Vector3>(degree, newCPs, newKnots);
            }

            return result;
        }

        public static NurbsCurve<Vector3> ToUnclampCurve(NurbsCurve<Vector3> curve)
        {
            int degree = curve.Degree;
            var knotVector = new List<double>(curve.Knots);
            var controlPoints = curve.ControlPoints;

            int n = controlPoints.Count - 1;
            var cw = new Vector4[controlPoints.Count];

            // Convert to homogeneous coordinates
            for (int i = 0; i <= n; i++)
            {
                var cp = controlPoints[i];
                cw[i] = new Vector4(cp.Value * (float)cp.Weight, (float)cp.Weight);
            }

            for (int i = 0; i <= degree - 2; i++)
            {
                knotVector[degree - i - 1] =
                    knotVector[degree - i] - (knotVector[n - i + 1] - knotVector[n - i]);
                int k = degree - 1;
                for (int j = i; j >= 0; j--)
                {
                    double alpha =
                        (knotVector[degree] - knotVector[k])
                        / (knotVector[degree + j + 1] - knotVector[k]);
                    cw[j] = (cw[j] - (float)alpha * cw[j + 1]) / (float)(1.0 - alpha);
                    k = k - 1;
                }
            }

            knotVector[0] =
                knotVector[1] - (knotVector[n - degree + 2] - knotVector[n - degree + 1]);

            for (int i = 0; i <= degree - 2; i++)
            {
                knotVector[n + i + 2] =
                    knotVector[n + i + 1] + (knotVector[degree + i + 1] - knotVector[degree + i]);
                for (int j = i; j >= 0; j--)
                {
                    double alpha =
                        (knotVector[n + 1] - knotVector[n - j])
                        / (knotVector[n - j + i + 2] - knotVector[n - j]);
                    cw[n - j] = (cw[n - j] - (float)(1.0 - alpha) * cw[n - j - 1]) / (float)alpha;
                }
            }

            knotVector[n + degree + 1] =
                knotVector[n + degree] + (knotVector[2 * degree] - knotVector[2 * degree - 1]);

            var newControlPoints = new ControlPoint<Vector3>[controlPoints.Count];
            for (int i = 0; i <= n; i++)
            {
                if (MathUtils.IsZero(cw[i].W))
                {
                    newControlPoints[i] = new ControlPoint<Vector3>(Vector3.Zero, 0.0);
                }
                else
                {
                    newControlPoints[i] = new ControlPoint<Vector3>(
                        new Vector3(cw[i].X, cw[i].Y, cw[i].Z) / cw[i].W,
                        cw[i].W
                    );
                }
            }

            return new NurbsCurve<Vector3>(degree, newControlPoints, knotVector);
        }

        public static double ApproximateLength(NurbsCurve<Vector3> curve, IntegratorType type)
        {
            if (IsLinear(curve))
            {
                var p0 = curve.ControlPoints[0].Value;
                var p1 = curve.ControlPoints[curve.ControlPoints.Count - 1].Value;
                return (p1 - p0).Length();
            }

            if (curve.Degree == 1)
            {
                double length = 0;
                for (int i = 0; i < curve.ControlPoints.Count - 1; i++)
                {
                    length += (
                        curve.ControlPoints[i + 1].Value - curve.ControlPoints[i].Value
                    ).Length();
                }
                return length;
            }

            var reCurve = Reparametrize(curve, 0.0, 1.0);
            int degree = reCurve.Degree;
            var knots = reCurve.Knots;
            var controlPoints = reCurve.ControlPoints;

            double totalLength = 0.0;
            switch (type)
            {
                case IntegratorType.Simpson:
                {
                    double start = knots[0];
                    double end = knots[knots.Count - 1];
                    totalLength = Integrator.Simpson(
                        (t, data) =>
                        {
                            var c = (NurbsCurve<Vector3>)data;
                            return ComputeRationalCurveDerivatives(c, 1, t)[1].Length();
                        },
                        reCurve,
                        start,
                        end
                    );
                    break;
                }
                case IntegratorType.GaussLegendre:
                {
                    var abscissae = Integrator.GaussLegendreAbscissae;
                    var weights = Integrator.GaussLegendreWeights;
                    int size = abscissae.Count;

                    for (int i = degree; i < knots.Count - degree - 1; i++)
                    {
                        double a = knots[i];
                        double b = knots[i + 1];
                        if (MathUtils.IsAlmostEqualTo(a, b))
                            continue;

                        double halfLen = (b - a) / 2.0;
                        double mid = (a + b) / 2.0;
                        double spanLength = 0.0;

                        for (int j = 0; j < size; j++)
                        {
                            double t = halfLen * abscissae[j] + mid;
                            var ders = ComputeRationalCurveDerivatives(reCurve, 1, t);
                            double derLen = ders[1].Length();
                            if (double.IsNaN(derLen))
                                derLen = 0.0;
                            spanLength += weights[j] * derLen;
                        }
                        totalLength += halfLen * spanLength;
                    }
                    break;
                }
                case IntegratorType.Chebyshev:
                {
                    var series = Integrator.ChebyshevSeries(32);
                    for (int i = degree; i < knots.Count - degree - 1; i++)
                    {
                        double a = knots[i];
                        double b = knots[i + 1];
                        if (MathUtils.IsAlmostEqualTo(a, b))
                            continue;

                        totalLength += Integrator.ClenshawCurtisQuadrature(
                            (t, data) =>
                            {
                                var c = (NurbsCurve<Vector3>)data;
                                return ComputeRationalCurveDerivatives(c, 1, t)[1].Length();
                            },
                            reCurve,
                            a,
                            b,
                            series,
                            Constants.DistanceEpsilon
                        );
                    }
                    break;
                }
            }
            return totalLength;
        }

        public static double GetParamOnCurve(
            NurbsCurve<Vector3> curve,
            double givenLength,
            IntegratorType type
        )
        {
            var knots = curve.Knots;
            double start = knots[0];
            double end = knots[knots.Count - 1];

            double totalLength = ApproximateLength(curve, type);
            if (MathUtils.IsLessThan(totalLength, givenLength, Constants.DistanceEpsilon))
            {
                return end;
            }
            if (MathUtils.IsAlmostEqualTo(givenLength, 0))
            {
                return start;
            }

            for (int i = 0; i < knots.Count; i++)
            {
                double knot = knots[i];
                if (MathUtils.IsAlmostEqualTo(knot, start) || MathUtils.IsAlmostEqualTo(knot, end))
                    continue;

                if (SplitAt(curve, knot, out var left, out var right))
                {
                    double length = ApproximateLength(left, type);
                    if (MathUtils.IsAlmostEqualTo(length, givenLength, Constants.DistanceEpsilon))
                    {
                        return knot;
                    }
                    if (MathUtils.IsGreaterThan(length, givenLength, Constants.DistanceEpsilon))
                    {
                        end = knot;
                        break;
                    }
                }
            }
            return GetParamByLength(curve, start, end, givenLength, type);
        }

        private static double GetParamByLength(
            NurbsCurve<Vector3> curve,
            double start,
            double end,
            double givenLength,
            IntegratorType type
        )
        {
            double low = start;
            double high = end;
            double mid = (low + high) / 2.0;

            for (int i = 0; i < 100; i++)
            {
                mid = (low + high) / 2.0;

                double length = Integrator.Simpson(
                    (t, data) =>
                    {
                        var c = (NurbsCurve<Vector3>)data;
                        return ComputeRationalCurveDerivatives(c, 1, t)[1].Length();
                    },
                    curve,
                    start,
                    mid
                );

                if (MathUtils.IsAlmostEqualTo(length, givenLength, Constants.DistanceEpsilon))
                    return mid;

                if (length < givenLength)
                    low = mid;
                else
                    high = mid;
            }
            return mid;
        }

        public static List<double> GetParamsOnCurve(
            NurbsCurve<Vector3> curve,
            double givenLength,
            IntegratorType type
        )
        {
            var result = new List<double>();

            var knots = curve.Knots;
            double end = knots[knots.Count - 1];

            double param = GetParamOnCurve(curve, givenLength, type);

            while (!MathUtils.IsAlmostEqualTo(param, end))
            {
                result.Add(param);

                if (!SplitAt(curve, param, out _, out var right))
                {
                    break;
                }

                param = GetParamOnCurve(right, givenLength, type);
            }
            return result;
        }

        public static bool IsClamp(NurbsCurve<Vector3> curve)
        {
            return KnotsUtils.IsClamped(curve.Degree, curve.Knots);
        }

        public static bool IsPeriodic(NurbsCurve<Vector3> curve)
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;
            int size = controlPoints.Count;

            if (KnotsUtils.IsClamped(degree, knots))
                return false;

            if (!KnotsUtils.IsUniform(knots))
                return false;

            if (size >= degree + degree)
            {
                bool flag = true;
                for (int i = 0; i < degree; i++)
                {
                    if (
                        !MathUtils.IsAlmostEqualTo(
                            controlPoints[i].Value,
                            controlPoints[size - degree + i].Value
                        )
                    )
                    {
                        flag = false;
                        break;
                    }
                }
                if (flag)
                {
                    return true;
                }
            }

            if (!IsClosed(curve))
                return false;

            double first = knots[0];
            double end = knots[knots.Count - 1];

            int cFirst = KnotsUtils.GetContinuity(degree, knots, first);
            int cEnd = KnotsUtils.GetContinuity(degree, knots, end);

            if (cFirst != cEnd)
                return false;

            var fDers = ComputeRationalCurveDerivatives(curve, cFirst, first);
            var eDers = ComputeRationalCurveDerivatives(curve, cEnd, end);

            for (int i = 0; i <= cFirst; i++)
            {
                var currentF = fDers[i];
                var currentE = eDers[i];

                var nf = Vector3.Normalize(currentF);
                var ne = Vector3.Normalize(currentE);

                bool hasSameDirection = MathUtils.IsAlmostEqualTo(nf, ne);
                bool hasSameMagnitude = MathUtils.IsAlmostEqualTo(
                    currentF.Length(),
                    currentE.Length()
                );
                if (!hasSameDirection || !hasSameMagnitude)
                {
                    return false;
                }
            }
            return true;
        }

        public static bool IsArc(NurbsCurve<Vector3> curve, out Vector3 center, out double radius)
        {
            center = Vector3.Zero;
            radius = 0.0;

            if (IsLinear(curve))
            {
                return false;
            }

            var knots = curve.Knots;
            double first = knots[0];
            double end = knots[knots.Count - 1];

            var p0 = GetPointOnCurve(curve, first);
            double param = IsClosed(curve) ? 0.5 * first + 0.5 * end : end;
            var p1 = GetPointOnCurve(curve, 0.5 * first + 0.5 * param);
            var p2 = GetPointOnCurve(curve, param);

            var v1 = p1 - p0;
            var v2 = p2 - p0;
            double v1v1 = Vector3.Dot(v1, v1);
            double v2v2 = Vector3.Dot(v2, v2);
            double v1v2 = Vector3.Dot(v1, v2);

            double det = v1v1 * v2v2 - v1v2 * v1v2;

            if (MathUtils.IsZero(det))
            {
                return false;
            }

            double baseVal = 0.5 / det;
            double k1 = baseVal * v2v2 * (v1v1 - v1v2);
            double k2 = baseVal * v1v1 * (v2v2 - v1v2);
            center = p0 + v1 * (float)k1 + v2 * (float)k2;
            radius = Vector3.Distance(center, p0);

            var (tessellatedPoints, _) = EquallyTessellate(curve);
            foreach (var point in tessellatedPoints)
            {
                double d = Vector3.Distance(point, center);
                if (!MathUtils.IsAlmostEqualTo(d, radius, Constants.DistanceEpsilon))
                {
                    return false;
                }
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

        public static void SplitArc(
            Vector3 start,
            Vector3 projectPoint,
            double projectPointWeight,
            Vector3 end,
            out Vector3 insertPointAtStartSide,
            out Vector3 splitPoint,
            out Vector3 insertPointAtEndSide,
            out double insertWeight
        )
        {
            insertPointAtStartSide = start + projectPoint;
            insertPointAtEndSide = end + projectPoint;
            splitPoint = (insertPointAtStartSide + insertPointAtEndSide) * 0.5f;
            insertWeight = Math.Sqrt(1 + projectPointWeight) * 0.5;
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

        public static int InsertKnot(
            NurbsCurve<Vector3> curve,
            double insertKnot,
            int times,
            out NurbsCurve<Vector3> result
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

            var updatedControlPoints = new ControlPoint<Vector3>[controlPoints.Count + times];
            for (int i = 0; i <= knotSpanIndex - degree; i++)
            {
                updatedControlPoints[i] = controlPoints[i];
            }
            for (int i = knotSpanIndex - originMultiplicity; i < controlPoints.Count; i++)
            {
                updatedControlPoints[i + times] = controlPoints[i];
            }

            var temp = new Vector4[degree - originMultiplicity + 1];
            for (int i = 0; i <= degree - originMultiplicity; i++)
            {
                var cp = controlPoints[knotSpanIndex - degree + i];
                temp[i] = new Vector4(cp.Value * (float)cp.Weight, (float)cp.Weight);
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
                updatedControlPoints[L] = MathUtils.IsZero(tVal.W)
                    ? new ControlPoint<Vector3>(Vector3.Zero, 0)
                    : new ControlPoint<Vector3>(
                        new Vector3(tVal.X, tVal.Y, tVal.Z) / tVal.W,
                        tVal.W
                    );

                if (degree - j - originMultiplicity > 0)
                {
                    var tBack = temp[degree - j - originMultiplicity];
                    int idx = knotSpanIndex + times - j - originMultiplicity;
                    updatedControlPoints[idx] = MathUtils.IsZero(tBack.W)
                        ? new ControlPoint<Vector3>(Vector3.Zero, 0)
                        : new ControlPoint<Vector3>(
                            new Vector3(tBack.X, tBack.Y, tBack.Z) / tBack.W,
                            tBack.W
                        );
                }
            }

            for (int i = L + 1; i < knotSpanIndex - originMultiplicity; i++)
            {
                var tVal = temp[i - L];
                updatedControlPoints[i] = MathUtils.IsZero(tVal.W)
                    ? new ControlPoint<Vector3>(Vector3.Zero, 0)
                    : new ControlPoint<Vector3>(
                        new Vector3(tVal.X, tVal.Y, tVal.Z) / tVal.W,
                        tVal.W
                    );
            }

            result = new NurbsCurve<Vector3>(
                degree,
                updatedControlPoints,
                insertedKnotVector.ToArray()
            );
            return times;
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

        public static bool Merge(
            NurbsCurve<Vector3> left,
            NurbsCurve<Vector3> right,
            out NurbsCurve<Vector3> result
        )
        {
            result = default;
            var cpL = left.ControlPoints;
            var cpR = right.ControlPoints;

            if (cpL.Count == 0 || cpR.Count == 0)
                return false;

            if (!MathUtils.IsAlmostEqualTo(cpL[cpL.Count - 1].Value, cpR[0].Value))
            {
                return false;
            }

            int degree = Math.Max(left.Degree, right.Degree);

            var tempL = Reparametrize(left, 0.0, 1.0);
            if (degree > left.Degree)
            {
                tempL = ElevateDegree(tempL, degree - left.Degree);
                tempL = Reparametrize(tempL, 0.0, 1.0);
            }

            var tempR = Reparametrize(right, 0.0, 1.0);
            if (degree > right.Degree)
            {
                tempR = ElevateDegree(tempR, degree - right.Degree);
                tempR = Reparametrize(tempR, 0.0, 1.0);
            }

            var kL = tempL.Knots;
            int lMulti = Polynomials.GetKnotMultiplicity(kL, kL[kL.Count - 1]);
            var kR = tempR.Knots;
            int rMulti = Polynomials.GetKnotMultiplicity(kR, kR[0]);

            if (lMulti != degree + 1 || rMulti != degree + 1)
            {
                return false;
            }

            var mergedPoints = new List<ControlPoint<Vector3>>(
                tempL.ControlPoints.Count + tempR.ControlPoints.Count
            );
            mergedPoints.AddRange(tempL.ControlPoints);
            mergedPoints.AddRange(tempR.ControlPoints);

            var shiftedR = Reparametrize(tempR, 1.0, 2.0);
            var kRShift = shiftedR.Knots;

            var mergedKnots = new List<double>(kL);
            for (int i = degree + 1; i < kRShift.Count; i++)
            {
                mergedKnots.Add(kRShift[i]);
            }

            var finalKnots = KnotsUtils.Rescale(mergedKnots, 0.0, 1.0);

            result = new NurbsCurve<Vector3>(degree, mergedPoints, finalKnots);
            return true;
        }
    }
}
