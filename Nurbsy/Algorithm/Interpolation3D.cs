using Stride.Core.Mathematics;

namespace Nurbsy.Algorithm
{
    internal static class Interpolation3D
    {
        public static double GetTotalChordLength(IReadOnlyList<Vector3> throughPoints)
        {
            double length = 0.0;
            int n = throughPoints.Count - 1;
            for (int i = 1; i <= n; i++)
            {
                length += Vector3.Distance(throughPoints[i], throughPoints[i - 1]);
            }
            return length;
        }

        public static double[] GetChordParameterization(IReadOnlyList<Vector3> throughPoints)
        {
            int size = throughPoints.Count;
            var uk = new double[size];
            uk[size - 1] = 1.0;

            double d = GetTotalChordLength(throughPoints);
            for (int i = 1; i < size - 1; i++)
            {
                uk[i] = uk[i - 1] + Vector3.Distance(throughPoints[i], throughPoints[i - 1]) / d;
            }
            return uk;
        }

        public static double GetCentripetalLength(IReadOnlyList<Vector3> throughPoints)
        {
            double length = 0.0;
            int n = throughPoints.Count - 1;
            for (int i = 1; i <= n; i++)
            {
                length += Math.Sqrt(Vector3.Distance(throughPoints[i], throughPoints[i - 1]));
            }
            return length;
        }

        public static double[] GetCentripetalParameterization(IReadOnlyList<Vector3> throughPoints)
        {
            int size = throughPoints.Count;
            var uk = new double[size];
            uk[size - 1] = 1.0;

            double d = GetCentripetalLength(throughPoints);
            for (int i = 1; i < size - 1; i++)
            {
                double dist = Vector3.Distance(throughPoints[i], throughPoints[i - 1]);
                uk[i] = uk[i - 1] + Math.Sqrt(dist) / d;
            }
            return uk;
        }

        public static bool ComputerWeightForRationalQuadraticInterpolation(
            Vector3 startPoint,
            Vector3 middleControlPoint,
            Vector3 endPoint,
            out double weight
        )
        {
            weight = 0.0;
            Vector3 SM = middleControlPoint - startPoint;
            Vector3 EM = middleControlPoint - endPoint;
            Vector3 SE = endPoint - startPoint;

            if (
                SM.LengthSquared() > MathUtil.ZeroTolerance
                && EM.LengthSquared() > MathUtil.ZeroTolerance
            )
            {
                Vector3 smNorm = Vector3.Normalize(SM);
                Vector3 emNorm = Vector3.Normalize(EM);

                if (
                    MathUtils.IsAlmostEqualTo(smNorm, emNorm)
                    || MathUtils.IsAlmostEqualTo(smNorm, -emNorm)
                )
                {
                    weight = 1.0;
                    return true;
                }
            }

            if (MathUtils.IsAlmostEqualTo(SM.Length(), EM.Length()))
            {
                double d = SM.Length() * SE.Length();
                if (MathUtils.IsZero(d))
                    return false;
                weight = Vector3.Dot(SM, SE) / d;
                return true;
            }

            Vector3 M = 0.5f * (startPoint + endPoint);
            Vector3 MR = middleControlPoint - M;
            double seLen = SE.Length();
            if (MathUtils.IsZero(seLen))
                return false;

            double ratio = SM.Length() / seLen;
            double frac = 1.0 / (ratio + 1.0);
            Vector3 D = endPoint + (float)frac * EM;
            Vector3 SD = D - startPoint;

            Vector3 dirSD =
                SD.LengthSquared() > MathUtil.ZeroTolerance ? Vector3.Normalize(SD) : Vector3.Zero;
            Vector3 dirMR =
                MR.LengthSquared() > MathUtil.ZeroTolerance ? Vector3.Normalize(MR) : Vector3.Zero;

            var type = Intersection3D.ComputeRays(
                startPoint,
                dirSD,
                M,
                dirMR,
                out double alf1,
                out double alf2,
                out Vector3 S1
            );
            if (type != CurveCurveIntersectionType.Intersecting)
                return false;

            double emLen = EM.Length();
            ratio = emLen / seLen;
            frac = 1.0 / (ratio + 1.0);
            Vector3 D2 = startPoint + (float)frac * SM;
            Vector3 ED = D2 - endPoint;

            Vector3 dirED =
                ED.LengthSquared() > MathUtil.ZeroTolerance ? Vector3.Normalize(ED) : Vector3.Zero;
            type = Intersection3D.ComputeRays(
                endPoint,
                dirED,
                M,
                dirMR,
                out alf1,
                out alf2,
                out Vector3 S2
            );
            if (type != CurveCurveIntersectionType.Intersecting)
                return false;

            Vector3 S = 0.5f * (S1 + S2);
            Vector3 MS = S - M;
            double msLen = MS.Length();
            double mrLen = MR.Length();
            if (MathUtils.IsZero(mrLen))
                return false;

            double s = msLen / mrLen;
            weight = s / (1.0 - s);
            return true;
        }

        public static bool GetSurfaceMeshParameterization(
            IReadOnlyList<IReadOnlyList<Vector3>> throughPoints,
            out double[] paramsU,
            out double[] paramsV
        )
        {
            int n = throughPoints.Count;
            if (n == 0)
            {
                paramsU = Array.Empty<double>();
                paramsV = Array.Empty<double>();
                return false;
            }
            int m = throughPoints[0].Count;

            var cds = new double[Math.Max(n, m)];
            paramsU = new double[n];
            paramsV = new double[m];

            int num = m;
            for (int l = 0; l < m; l++)
            {
                double total = 0.0;
                for (int k = 1; k < n; k++)
                {
                    double dist = Vector3.Distance(throughPoints[k][l], throughPoints[k - 1][l]);
                    cds[k] = dist;
                    total += dist;
                }

                if (MathUtils.IsAlmostEqualTo(total, 0.0))
                {
                    num--;
                }
                else
                {
                    double d = 0.0;
                    for (int k = 1; k < n; k++)
                    {
                        d += cds[k];
                        paramsU[k] += d / total;
                    }
                }
            }

            if (num == 0)
                return false;

            for (int k = 1; k < n - 1; k++)
            {
                paramsU[k] /= num;
            }
            paramsU[n - 1] = 1.0;

            num = n;
            for (int k = 0; k < n; k++)
            {
                double total = 0.0;
                for (int l = 1; l < m; l++)
                {
                    double dist = Vector3.Distance(throughPoints[k][l], throughPoints[k][l - 1]);
                    cds[l] = dist;
                    total += dist;
                }

                if (MathUtils.IsAlmostEqualTo(total, 0.0))
                {
                    num--;
                }
                else
                {
                    double d = 0.0;
                    for (int l = 1; l < m; l++)
                    {
                        d += cds[l];
                        paramsV[l] += d / total;
                    }
                }
            }

            if (num == 0)
                return false;

            for (int l = 1; l < m - 1; l++)
            {
                paramsV[l] /= num;
            }
            paramsV[m - 1] = 1.0;

            return true;
        }

        public static Vector3[] ComputeTangents(IReadOnlyList<Vector3> throughPoints)
        {
            int size = throughPoints.Count;
            var tangents = new Vector3[size];
            var qq = new Vector3[size];
            var delta = new double[size];

            var paramsT = GetChordParameterization(throughPoints);

            for (int i = 1; i < size; i++)
            {
                delta[i] = paramsT[i] - paramsT[i - 1];
                qq[i] = throughPoints[i] - throughPoints[i - 1];
            }

            for (int i = 1; i < size - 1; i++)
            {
                double den = delta[i] + delta[i + 1];
                if (MathUtils.IsZero(den))
                {
                    tangents[i] = Vector3.Normalize(qq[i + 1]);
                    continue;
                }
                double a = delta[i] / den;
                Vector3 vec = (float)(1.0 - a) * qq[i] + (float)a * qq[i + 1];
                tangents[i] = Vector3.Normalize(vec);
            }

            if (size > 1)
            {
                if (!MathUtils.IsZero(delta[1]))
                    tangents[0] = Vector3.Normalize(2.0f * qq[1] / (float)delta[1] - tangents[1]);
                else
                    tangents[0] = tangents[1];

                if (!MathUtils.IsZero(delta[size - 1]))
                {
                    tangents[size - 1] = Vector3.Normalize(
                        2.0f * qq[size - 1] / (float)delta[size - 1] - tangents[size - 2]
                    );
                }
                else
                    tangents[size - 1] = tangents[size - 2];
            }
            else if (size == 1)
            {
                tangents[0] = Vector3.UnitX;
            }

            return tangents;
        }

        public static bool TryComputeTangents(
            IReadOnlyList<Vector3> throughPoints,
            out Vector3[] tangents
        )
        {
            int size = throughPoints.Count;
            if (size < 5)
            {
                tangents = null;
                return false;
            }

            tangents = new Vector3[size];
            for (int k = 2; k < size - 2; k++)
            {
                Vector3 qk_1 = Getqk(throughPoints, k - 1);
                Vector3 qk = Getqk(throughPoints, k);
                Vector3 qk1 = Getqk(throughPoints, k + 1);
                Vector3 qk2 = Getqk(throughPoints, k + 2);
                tangents[k] = GetTk(qk_1, qk, qk1, qk2);
            }

            int n = size - 1;
            Vector3 q0 = 2.0f * Getqk(throughPoints, 1) - Getqk(throughPoints, 2);
            Vector3 q_1 = 2.0f * q0 - Getqk(throughPoints, 1);
            Vector3 qn1 = 2.0f * Getqk(throughPoints, n) - Getqk(throughPoints, n - 1);
            Vector3 qn2 = 2.0f * qn1 - Getqk(throughPoints, n);

            tangents[0] = GetTk(q_1, q0, Getqk(throughPoints, 1), Getqk(throughPoints, 2));
            tangents[1] = GetTk(
                q0,
                Getqk(throughPoints, 1),
                Getqk(throughPoints, 2),
                Getqk(throughPoints, 3)
            );
            tangents[size - 2] = GetTk(
                Getqk(throughPoints, n - 2),
                Getqk(throughPoints, n - 1),
                Getqk(throughPoints, n),
                qn1
            );
            tangents[size - 1] = GetTk(
                Getqk(throughPoints, n - 1),
                Getqk(throughPoints, n),
                qn1,
                qn2
            );

            return true;
        }

        private static Vector3 Getqk(IReadOnlyList<Vector3> pts, int index)
        {
            return pts[index] - pts[index - 1];
        }

        private static double Getak(Vector3 qk_1, Vector3 qk, Vector3 qk1, Vector3 qk2)
        {
            double cp1 = Vector3.Cross(qk_1, qk).Length();
            double cp2 = Vector3.Cross(qk1, qk2).Length();
            double sum = cp1 + cp2;
            if (MathUtils.IsZero(sum))
                return 0.5;
            return cp1 / sum;
        }

        private static Vector3 GetTk(Vector3 qk_1, Vector3 qk, Vector3 qk1, Vector3 qk2)
        {
            double ak = Getak(qk_1, qk, qk1, qk2);
            Vector3 v = (float)(1.0 - ak) * qk + (float)ak * qk1;
            return Vector3.Normalize(v);
        }
    }
}
