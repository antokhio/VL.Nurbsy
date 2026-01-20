using Stride.Core.Mathematics;

namespace Nurbsy.Algorithm
{
    /// <summary>
    /// Helper methods for computing surface area using numerical integration.
    /// </summary>
    public static class SurfaceAreaHelper
    {
        /// <summary>
        /// Compute the area element (first fundamental form determinant) at a UV point.
        /// ds = sqrt(E*G - F*F) where E, F, G are first fundamental form coefficients.
        /// </summary>
        public static double ComputeAreaElement(Vector3 Su, Vector3 Sv)
        {
            double E = Vector3.Dot(Su, Su);
            double F = Vector3.Dot(Su, Sv);
            double G = Vector3.Dot(Sv, Sv);
            double det = E * G - F * F;
            return det > 0 ? Math.Sqrt(det) : 0.0;
        }

        /// <summary>
        /// Compute the area element for 2D surfaces.
        /// </summary>
        public static double ComputeAreaElement(Vector2 Su, Vector2 Sv)
        {
            double E = Vector2.Dot(Su, Su);
            double F = Vector2.Dot(Su, Sv);
            double G = Vector2.Dot(Sv, Sv);
            double det = E * G - F * F;
            return det > 0 ? Math.Sqrt(det) : 0.0;
        }

        /// <summary>
        /// Compute approximate area using Simpson's adaptive rule.
        /// </summary>
        public static double ApproximateAreaSimpson<T>(
            in NurbsSurface<T> surface,
            Func<NurbsSurface<T>, double, double, double> areaElementFunc,
            double tolerance = 1e-6
        )
            where T : struct
        {
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;

            double startU = knotsU[0];
            double endU = knotsU[^1];
            double startV = knotsV[0];
            double endV = knotsV[^1];

            // Reparametrize to [0,1] x [0,1] for consistent integration
            double scaleU = endU - startU;
            double scaleV = endV - startV;

            double area = 0.0;

            // Adaptive subdivision using a stack
            var stack = new Stack<(double u1, double u2, double v1, double v2, double estArea)>();

            // Initial estimate over entire domain
            double initArea = SimpsonIntegrate(
                surface,
                areaElementFunc,
                startU,
                endU,
                startV,
                endV
            );
            stack.Push((startU, endU, startV, endV, initArea));

            while (stack.Count > 0)
            {
                var (u1, u2, v1, v2, estArea) = stack.Pop();

                double du = u2 - u1;
                double dv = v2 - v1;
                double hdu = 0.5 * du;
                double hdv = 0.5 * dv;

                // Bisect into 4 parts
                double area1 = SimpsonIntegrate(
                    surface,
                    areaElementFunc,
                    u1,
                    u1 + hdu,
                    v1,
                    v1 + hdv
                );
                double area2 = SimpsonIntegrate(
                    surface,
                    areaElementFunc,
                    u1 + hdu,
                    u2,
                    v1,
                    v1 + hdv
                );
                double area3 = SimpsonIntegrate(
                    surface,
                    areaElementFunc,
                    u1,
                    u1 + hdu,
                    v1 + hdv,
                    v2
                );
                double area4 = SimpsonIntegrate(
                    surface,
                    areaElementFunc,
                    u1 + hdu,
                    u2,
                    v1 + hdv,
                    v2
                );

                double areaNew = area1 + area2 + area3 + area4;

                if (Math.Abs(areaNew - estArea) < tolerance)
                {
                    area += areaNew;
                }
                else
                {
                    stack.Push((u1, u1 + hdu, v1, v1 + hdv, area1));
                    stack.Push((u1 + hdu, u2, v1, v1 + hdv, area2));
                    stack.Push((u1, u1 + hdu, v1 + hdv, v2, area3));
                    stack.Push((u1 + hdu, u2, v1 + hdv, v2, area4));
                }
            }

            return area;
        }

        /// <summary>
        /// Simple Simpson's rule integration over a rectangular domain.
        /// </summary>
        private static double SimpsonIntegrate<T>(
            in NurbsSurface<T> surface,
            Func<NurbsSurface<T>, double, double, double> areaElementFunc,
            double u1,
            double u2,
            double v1,
            double v2
        )
            where T : struct
        {
            double du = u2 - u1;
            double dv = v2 - v1;
            double um = (u1 + u2) * 0.5;
            double vm = (v1 + v2) * 0.5;

            // 9-point Simpson's rule
            double f00 = areaElementFunc(surface, u1, v1);
            double f01 = areaElementFunc(surface, u1, vm);
            double f02 = areaElementFunc(surface, u1, v2);
            double f10 = areaElementFunc(surface, um, v1);
            double f11 = areaElementFunc(surface, um, vm);
            double f12 = areaElementFunc(surface, um, v2);
            double f20 = areaElementFunc(surface, u2, v1);
            double f21 = areaElementFunc(surface, u2, vm);
            double f22 = areaElementFunc(surface, u2, v2);

            return (du * dv / 36.0)
                * (f00 + f22 + f02 + f20 + 4.0 * (f01 + f10 + f12 + f21) + 16.0 * f11);
        }

        /// <summary>
        /// Compute approximate area using Gauss-Legendre quadrature.
        /// Decomposes surface to Bezier patches and integrates each.
        /// </summary>
        public static double ApproximateAreaGaussLegendre<T>(
            in NurbsSurface<T> surface,
            Func<NurbsSurface<T>, double, double, double> areaElementFunc
        )
            where T : struct
        {
            double area = 0.0;

            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;

            // Get unique knot spans
            var uniqueKnotsU = KnotsUtils.GetUniqueKnots(knotsU);
            var uniqueKnotsV = KnotsUtils.GetUniqueKnots(knotsV);

            // Gauss-Legendre abscissae and weights (5-point)
            var abscissae = Integrator.GaussLegendreAbscissae;
            var weights = Integrator.GaussLegendreWeights;
            int n = abscissae.Count;

            // Integrate over each knot span
            for (int i = 0; i < uniqueKnotsU.Count - 1; i++)
            {
                double a = uniqueKnotsU[i];
                double b = uniqueKnotsU[i + 1];
                double coeffU = (b - a) * 0.5;
                double midU = (a + b) * 0.5;

                for (int j = 0; j < uniqueKnotsV.Count - 1; j++)
                {
                    double c = uniqueKnotsV[j];
                    double d = uniqueKnotsV[j + 1];
                    double coeffV = (d - c) * 0.5;
                    double midV = (c + d) * 0.5;

                    double patchArea = 0.0;

                    for (int iu = 0; iu < n; iu++)
                    {
                        double u = coeffU * abscissae[iu] + midU;

                        for (int iv = 0; iv < n; iv++)
                        {
                            double v = coeffV * abscissae[iv] + midV;
                            double ds = areaElementFunc(surface, u, v);
                            patchArea += weights[iu] * weights[iv] * ds;
                        }
                    }

                    patchArea *= coeffU * coeffV;
                    area += patchArea;
                }
            }

            return area;
        }

        /// <summary>
        /// Compute approximate area using Chebyshev/Clenshaw-Curtis quadrature.
        /// </summary>
        public static double ApproximateAreaChebyshev<T>(
            in NurbsSurface<T> surface,
            Func<NurbsSurface<T>, double, double, double> areaElementFunc,
            int order = 16
        )
            where T : struct
        {
            double area = 0.0;

            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;

            // Get unique knot spans
            var uniqueKnotsU = KnotsUtils.GetUniqueKnots(knotsU);
            var uniqueKnotsV = KnotsUtils.GetUniqueKnots(knotsV);

            // Chebyshev nodes and weights
            var (nodes, weights) = GetClenshawCurtisNodesAndWeights(order);

            // Integrate over each knot span
            for (int i = 0; i < uniqueKnotsU.Count - 1; i++)
            {
                double a = uniqueKnotsU[i];
                double b = uniqueKnotsU[i + 1];
                double coeffU = (b - a) * 0.5;
                double midU = (a + b) * 0.5;

                for (int j = 0; j < uniqueKnotsV.Count - 1; j++)
                {
                    double c = uniqueKnotsV[j];
                    double d = uniqueKnotsV[j + 1];
                    double coeffV = (d - c) * 0.5;
                    double midV = (c + d) * 0.5;

                    double patchArea = 0.0;

                    for (int iu = 0; iu < order; iu++)
                    {
                        double u = coeffU * nodes[iu] + midU;

                        for (int iv = 0; iv < order; iv++)
                        {
                            double v = coeffV * nodes[iv] + midV;
                            double ds = areaElementFunc(surface, u, v);
                            patchArea += weights[iu] * weights[iv] * ds;
                        }
                    }

                    patchArea *= coeffU * coeffV;
                    area += patchArea;
                }
            }

            return area;
        }

        private static (double[] nodes, double[] weights) GetClenshawCurtisNodesAndWeights(int n)
        {
            var nodes = new double[n];
            var weights = new double[n];

            for (int k = 0; k < n; k++)
            {
                nodes[k] = Math.Cos(Math.PI * k / (n - 1));

                double w = 1.0;
                for (int j = 1; j <= (n - 1) / 2; j++)
                {
                    double b = j == (n - 1) / 2 && (n - 1) % 2 == 0 ? 1.0 : 2.0;
                    w -= b * Math.Cos(2.0 * j * Math.PI * k / (n - 1)) / (4.0 * j * j - 1.0);
                }

                if (k == 0 || k == n - 1)
                    w /= (n - 1);
                else
                    w *= 2.0 / (n - 1);

                weights[k] = w;
            }

            return (nodes, weights);
        }
    }
}
