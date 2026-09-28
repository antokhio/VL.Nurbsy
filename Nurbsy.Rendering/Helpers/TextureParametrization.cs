using Nurbsy.Algorithm;

namespace Nurbsy.Rendering.Helpers
{
    /// <summary>
    /// Maps normalized texture coordinates to surface parameters, one direction at a time.
    /// <para>
    /// The texture coordinate of a parameter <c>u</c> is defined as
    /// <c>tau(u) = sum_i N_i,p(u) * i / (n - 1)</c>, i.e. the B-spline whose coefficients are the
    /// normalized control point indices. The mapping depends only on the degree, the knot vector
    /// and the control point count, never on control point positions, so:
    /// </para>
    /// <list type="bullet">
    /// <item>moving a control point drags its region of the texture with it,</item>
    /// <item>an evenly spaced control grid yields an exactly uniform texture and vertex
    /// distribution (no compression at clamped ends for degree &gt; 1), because by linear
    /// precision the surface is then an affine function of <c>tau</c>.</item>
    /// </list>
    /// </summary>
    public static class TextureParametrization
    {
        private const int MaxIterations = 60;

        /// <summary>
        /// Resolves the surface parameters for <paramref name="samples"/> evenly spaced
        /// texture coordinates <c>t = k / (samples - 1)</c>.
        /// </summary>
        public static float[] Resolve(
            int degree,
            IReadOnlyList<double> knots,
            int controlPointCount,
            int samples
        )
        {
            samples = Math.Max(2, samples);
            var result = new float[samples];

            double min = knots[0];
            double max = knots[knots.Count - 1];
            double tauMin = Evaluate(degree, knots, controlPointCount, min);
            double tauMax = Evaluate(degree, knots, controlPointCount, max);
            double range = tauMax - tauMin;

            result[0] = (float)min;
            result[samples - 1] = (float)max;

            if (range <= 0)
            {
                for (int k = 1; k < samples - 1; k++)
                    result[k] = (float)(min + (max - min) * k / (samples - 1));
                return result;
            }

            double lo = min;
            for (int k = 1; k < samples - 1; k++)
            {
                double target = tauMin + range * k / (samples - 1);

                // tau is monotonically non-decreasing, so the previous solution bounds the search.
                double a = lo;
                double b = max;
                for (int it = 0; it < MaxIterations && b - a > 1e-12; it++)
                {
                    double mid = 0.5 * (a + b);
                    if (Evaluate(degree, knots, controlPointCount, mid) < target)
                        a = mid;
                    else
                        b = mid;
                }

                lo = 0.5 * (a + b);
                result[k] = (float)lo;
            }

            return result;
        }

        /// <summary>
        /// Evaluates the (unnormalized) texture coordinate <c>tau(u)</c>.
        /// </summary>
        public static double Evaluate(
            int degree,
            IReadOnlyList<double> knots,
            int controlPointCount,
            double u
        )
        {
            int span = Polynomials.GetKnotSpanIndex(degree, knots, u);
            double[] basis = Polynomials.BasisFunctions(span, degree, knots, u);

            double denominator = Math.Max(1, controlPointCount - 1);
            double tau = 0.0;
            for (int k = 0; k <= degree; k++)
            {
                tau += basis[k] * (span - degree + k) / denominator;
            }
            return tau;
        }
    }
}
