using Stride.Core.Mathematics;

namespace Nurbsy.Algorithm
{
    internal static class Interpolation3D
    {
        public static double GetTotalChordLength(IReadOnlyList<Vector3> throughPoints)
        {
            int n = throughPoints.Count - 1;
            double length = 0.0;
            for (int i = 1; i <= n; i++)
            {
                length += Vector3.Distance(throughPoints[i], throughPoints[i - 1]);
            }
            return length;
        }

        public static double[] GetChordParameterization(IReadOnlyList<Vector3> throughPoints)
        {
            int size = throughPoints.Count;
            int n = size - 1;

            var uk = new double[size];
            uk[0] = 0.0;
            uk[n] = 1.0;

            double d = GetTotalChordLength(throughPoints);
            if (MathUtils.IsZero(d))
            {
                for (int i = 1; i < n; i++)
                    uk[i] = (double)i / n;
            }
            else
            {
                for (int i = 1; i < n; i++)
                {
                    uk[i] =
                        uk[i - 1] + Vector3.Distance(throughPoints[i], throughPoints[i - 1]) / d;
                }
            }
            return uk;
        }

        public static double GetCentripetalLength(IReadOnlyList<Vector3> throughPoints)
        {
            int n = throughPoints.Count - 1;
            double length = 0.0;
            for (int i = 1; i <= n; i++)
            {
                length += Math.Sqrt(Vector3.Distance(throughPoints[i], throughPoints[i - 1]));
            }
            return length;
        }

        public static double[] GetCentripetalParameterization(IReadOnlyList<Vector3> throughPoints)
        {
            int size = throughPoints.Count;
            int n = size - 1;

            var uk = new double[size];
            uk[0] = 0.0;
            uk[n] = 1.0;

            double d = GetCentripetalLength(throughPoints);
            if (MathUtils.IsZero(d))
            {
                for (int i = 1; i < n; i++)
                    uk[i] = (double)i / n;
            }
            else
            {
                for (int i = 1; i < n; i++)
                {
                    uk[i] =
                        uk[i - 1]
                        + Math.Sqrt(Vector3.Distance(throughPoints[i], throughPoints[i - 1])) / d;
                }
            }
            return uk;
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
                var qk_1 = Getqk(throughPoints, k - 1);
                var qk = Getqk(throughPoints, k);
                var qk1 = Getqk(throughPoints, k + 1);
                var qk2 = Getqk(throughPoints, k + 2);

                tangents[k] = GetTk(qk_1, qk, qk1, qk2);
            }

            int n = size - 1;
            var q0 = 2.0f * Getqk(throughPoints, 1) - Getqk(throughPoints, 2);
            var q_1 = 2.0f * q0 - Getqk(throughPoints, 1);
            var qn1 = 2.0f * Getqk(throughPoints, n) - Getqk(throughPoints, n - 1);
            var qn2 = 2.0f * qn1 - Getqk(throughPoints, n);

            tangents[0] = GetTk(q_1, q0, Getqk(throughPoints, 1), Getqk(throughPoints, 2));
            tangents[1] = GetTk(
                q0,
                Getqk(throughPoints, 1),
                Getqk(throughPoints, 2),
                Getqk(throughPoints, 3)
            );
            tangents[n - 1] = GetTk(
                Getqk(throughPoints, n - 2),
                Getqk(throughPoints, n - 1),
                Getqk(throughPoints, n),
                qn1
            );
            tangents[n] = GetTk(Getqk(throughPoints, n - 1), Getqk(throughPoints, n), qn1, qn2);

            return true;
        }

        public static Vector3[] ComputeTangents(IReadOnlyList<Vector3> throughPoints)
        {
            int size = throughPoints.Count;
            var tangents = new Vector3[size];
            var delta = new double[size];
            var qq = new Vector3[size];

            var paramsT = GetChordParameterization(throughPoints);

            for (int i = 1; i < size; i++)
            {
                delta[i] = paramsT[i] - paramsT[i - 1];
                qq[i] = throughPoints[i] - throughPoints[i - 1];
            }

            for (int i = 1; i < size - 1; i++)
            {
                double a = delta[i] / (delta[i] + delta[i + 1]);
                tangents[i] = Vector3.Normalize((float)(1.0 - a) * qq[i] + (float)a * qq[i + 1]);
            }

            if (delta[1] > MathUtils.Epsilon)
                tangents[0] = Vector3.Normalize(2.0f * qq[1] / (float)delta[1] - tangents[1]);
            else
                tangents[0] = Vector3.Zero;

            if (delta[size - 1] > MathUtils.Epsilon)
                tangents[size - 1] = Vector3.Normalize(
                    2.0f * qq[size - 1] / (float)delta[size - 1] - tangents[size - 2]
                );
            else
                tangents[size - 1] = Vector3.Zero;

            return tangents;
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

            Vector3 SMNorm = Vector3.Normalize(SM);
            Vector3 EMNorm = Vector3.Normalize(EM);

            if (
                MathUtils.IsAlmostEqualTo(SMNorm, EMNorm)
                || MathUtils.IsAlmostEqualTo(SMNorm, -EMNorm)
            )
            {
                weight = 1.0;
                return true;
            }

            if (MathUtils.IsAlmostEqualTo(SM.Length(), EM.Length()))
            {
                weight = Vector3.Dot(SM, SE) / (SM.Length() * SE.Length());
                return true;
            }

            Vector3 M = 0.5f * (startPoint + endPoint);
            Vector3 MR = middleControlPoint - M;
            double ratio = SM.Length() / SE.Length();
            double frac = 1.0 / (ratio + 1.0);
            Vector3 D = endPoint + (float)frac * EM;
            Vector3 SD = D - startPoint;

            var type1 = Intersection3D.ComputeRays(
                startPoint,
                Vector3.Normalize(SD),
                M,
                Vector3.Normalize(MR),
                out double alf1,
                out double alf2,
                out Vector3 S1
            );
            if (type1 != CurveCurveIntersectionType.Intersecting)
                return false;

            ratio = EM.Length() / SE.Length();
            frac = 1.0 / (ratio + 1.0);
            D = startPoint + (float)frac * SM;
            Vector3 ED = D - endPoint;

            var type2 = Intersection3D.ComputeRays(
                endPoint,
                Vector3.Normalize(ED),
                M,
                Vector3.Normalize(MR),
                out alf1,
                out alf2,
                out Vector3 S2
            );
            if (type2 != CurveCurveIntersectionType.Intersecting)
                return false;

            Vector3 S = 0.5f * (S1 + S2);
            Vector3 MS = S - M;
            double s = MS.Length() / MR.Length();
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
            int m = throughPoints[0].Count;

            paramsU = new double[n];
            paramsV = new double[m];
            var cds = new double[Math.Max(n, m)];

            int num = m;
            for (int l = 0; l < m; l++)
            {
                double total = 0.0;
                for (int k = 1; k < n; k++)
                {
                    cds[k] = Vector3.Distance(throughPoints[k][l], throughPoints[k - 1][l]);
                    total += cds[k];
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
                    cds[l] = Vector3.Distance(throughPoints[k][l], throughPoints[k][l - 1]);
                    total += cds[l];
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

        private static Vector3 Getqk(IReadOnlyList<Vector3> throughPoints, int index)
        {
            return throughPoints[index] - throughPoints[index - 1];
        }

        private static double Getak(Vector3 qk_1, Vector3 qk, Vector3 qk1, Vector3 qk2)
        {
            double cross1 = Vector3.Cross(qk_1, qk).Length();
            double cross2 = Vector3.Cross(qk1, qk2).Length();
            double den = cross1 + cross2;
            if (MathUtils.IsZero(den))
                return 0.0;
            return cross1 / den;
        }

        private static Vector3 GetTk(Vector3 qk_1, Vector3 qk, Vector3 qk1, Vector3 qk2)
        {
            double ak = Getak(qk_1, qk, qk1, qk2);
            return Vector3.Normalize((float)(1.0 - ak) * qk + (float)ak * qk1);
        }
    }
}
