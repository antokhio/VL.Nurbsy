/*
 * Ported by:
 * 2025 - Anton Kalabuhov (antokhio)
 *
 * Original Author:
 * 2023 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 * Based on LNLib: https://github.com/BIMCoderLiang/LNLib
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

namespace Nurbsy.Algorithm
{
    /// <summary>
    /// Provides low-level algorithms for polynomial evaluation, B-Spline basis functions, and matrix conversions.
    /// </summary>
    public static class Polynomials
    {
        /// <summary>
        /// Evaluates a polynomial using Horner's method.
        /// </summary>
        /// <param name="degree">The degree of the polynomial.</param>
        /// <param name="coefficients">The coefficients of the polynomial.</param>
        /// <param name="paramT">The parameter value at which to evaluate.</param>
        /// <returns>The value of the polynomial at <paramref name="paramT"/>.</returns>
        public static double Horner(int degree, IReadOnlyList<double> coefficients, double paramT)
        {
            Validate.Argument(degree > 0, nameof(degree), "Degree must be greater than zero.");
            Validate.Argument(
                Validate.IsValidBezier(degree, coefficients.Count),
                nameof(degree),
                "Coefficients size equal degree plus one."
            );
            double result = coefficients[degree];
            for (int i = degree - 1; i >= 0; i--)
            {
                result = result * paramT + coefficients[i];
            }

            return result;
        }

        /// <summary>
        /// Evaluates a bivariate polynomial (surface) using Horner's method.
        /// </summary>
        /// <param name="degreeU">The degree in the U direction.</param>
        /// <param name="degreeV">The degree in the V direction.</param>
        /// <param name="coefficients">The grid of coefficients.</param>
        /// <param name="u">The U parameter.</param>
        /// <param name="v">The V parameter.</param>
        /// <returns>The evaluated value.</returns>
        public static double Horner(
            int degreeU,
            int degreeV,
            IReadOnlyList<IReadOnlyList<double>> coefficients,
            double u,
            double v
        )
        {
            Validate.Argument(degreeU <= 0, nameof(degreeU), "DegreeU must be greater than zero.");
            Validate.Argument(degreeV <= 0, nameof(degreeV), "DegreeV must be greater than zero.");
            Validate.Argument(
                coefficients.Count == 0,
                nameof(coefficients),
                "Coefficients size must be greater than zero."
            );
            Validate.Argument(
                !Validate.IsValidBezier(degreeU, coefficients.Count),
                nameof(degreeU),
                "Coefficients row size equals degreeU plus one."
            );
            Validate.Argument(
                !Validate.IsValidBezier(degreeV, coefficients[0].Count),
                nameof(degreeV),
                "Coefficients column size equals degreeV plus one."
            );
            Validate.Range(u, 0.0, 1.0, nameof(u));
            Validate.Range(v, 0.0, 1.0, nameof(v));

            var temp = new double[degreeU + 1];
            for (int i = 0; i <= degreeU; i++)
            {
                temp[i] = Horner(degreeV, coefficients[i], v);
            }
            return Horner(degreeU, temp, u);
        }

        /// <summary>
        /// Computes the value of the i-th Bernstein polynomial of degree p at t.
        /// </summary>
        /// <param name="index">The index i of the Bernstein polynomial.</param>
        /// <param name="degree">The degree p of the polynomial.</param>
        /// <param name="paramT">The parameter t.</param>
        /// <returns>The value B_{i,p}(t).</returns>
        public static double Bernstein(int index, int degree, double paramT)
        {
            Validate.Argument(
                index < 0 || index > degree,
                nameof(index),
                "Index must between zero and degree."
            );
            Validate.Argument(
                degree < 0,
                nameof(degree),
                "Degree must be greater than or equal zero."
            );
            Validate.Range(paramT, 0.0, 1.0, nameof(paramT));

            if (degree == 0)
                return 1.0;

            if (index == 0)
                return Math.Pow(1.0 - paramT, degree);

            if (index == degree)
                return Math.Pow(paramT, degree);

            var temp = new double[degree + 1];
            temp[degree - index] = 1.0;
            double t1 = 1.0 - paramT;

            for (int k = index; k <= degree; k++)
            {
                for (int j = degree; j >= k; j--)
                {
                    temp[j] = t1 * temp[j] + paramT * temp[j - 1];
                }
            }

            return temp[degree];
        }

        /// <summary>
        /// Computes values of all Bernstein polynomials of degree p at t.
        /// </summary>
        /// <param name="degree">The degree p.</param>
        /// <param name="paramT">The parameter t.</param>
        /// <returns>An array containing B_{0,p}(t), ..., B_{p,p}(t).</returns>
        public static double[] AllBernstein(int degree, double paramT)
        {
            Validate.Argument(
                degree < 0,
                nameof(degree),
                "Degree must be greater than zero or equal to zero."
            );
            Validate.Range(paramT, 0.0, 1.0, nameof(paramT));

            var bernsteinArray = new double[degree + 1];
            bernsteinArray[0] = 1.0;

            if (degree == 0)
            {
                return bernsteinArray;
            }

            double t1 = 1.0 - paramT;
            for (int j = 1; j <= degree; j++)
            {
                double saved = 0.0;
                for (int k = 0; k < j; k++)
                {
                    double temp = bernsteinArray[k];
                    bernsteinArray[k] = saved + t1 * temp;
                    saved = paramT * temp;
                }
                bernsteinArray[j] = saved;
            }

            return bernsteinArray;
        }

        /// <summary>
        /// Calculates the multiplicity of a knot in the knot vector.
        /// </summary>
        /// <param name="knots">The knot vector.</param>
        /// <param name="knot">The knot value.</param>
        /// <returns>The number of times the knot appears in the vector.</returns>
        public static int GetKnotMultiplicity(IReadOnlyList<double> knots, double knot)
        {
            Validate.Argument(
                knots.Count > 0,
                nameof(knots),
                "Knots size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knots),
                nameof(knots),
                "Knots must be a nondecreasing sequence of real numbers."
            );
            Validate.Range(knot, knots[0], knots[knots.Count - 1], nameof(knot));

            int size = knots.Count;
            int multi = 0;

            for (int index = 0; index < size; index++)
            {
                if (MathUtils.IsAlmostEqualTo(knot, knots[index]))
                {
                    multi++;
                }
            }

            return multi;
        }

        /// <summary>
        /// Finds the knot span index corresponding to a parameter value.
        /// </summary>
        /// <param name="degree">The degree of the basis functions.</param>
        /// <param name="knots">The knot vector.</param>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>The index of the knot span.</returns>
        public static int GetKnotSpanIndex(int degree, IReadOnlyList<double> knots, double paramT)
        {
            Validate.Argument(
                degree >= 0,
                nameof(degree),
                "Degree must be greater than or equal zero."
            );
            Validate.Argument(
                knots.Count > 0,
                nameof(knots),
                "Knots size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knots),
                nameof(knots),
                "Knots must be a nondecreasing sequence of real numbers."
            );
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            int n = knots.Count - degree - 2;
            if (MathUtils.IsGreaterThanOrEqual(paramT, knots[n + 1]))
            {
                return n;
            }
            if (MathUtils.IsLessThanOrEqual(paramT, knots[degree]))
            {
                return degree;
            }

            int low = degree;
            int high = n + 1;
            int mid = (low + high) / 2;

            while (paramT < knots[mid] || paramT >= knots[mid + 1])
            {
                if (paramT < knots[mid])
                {
                    high = mid;
                }
                else
                {
                    low = mid;
                }
                mid = (low + high) / 2;
            }
            return mid;
        }

        /// <summary>
        /// Computes the non-vanishing basis functions at a given parameter value.
        /// </summary>
        /// <param name="spanIndex">The knot span index.</param>
        /// <param name="degree">The degree of the basis functions.</param>
        /// <param name="knots">The knot vector.</param>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>Array of basis function values N_{i-p, p}(t) ... N_{i, p}(t).</returns>
        public static double[] BasisFunctions(
            int spanIndex,
            int degree,
            IReadOnlyList<double> knots,
            double paramT
        )
        {
            Validate.Argument(
                spanIndex >= 0,
                nameof(spanIndex),
                "SpanIndex must be greater than or equal zero."
            );
            Validate.Argument(
                degree >= 0,
                nameof(degree),
                "Degree must be greater than or equal zero."
            );
            Validate.Argument(
                knots.Count > 0,
                nameof(knots),
                "Knots size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knots),
                nameof(knots),
                "Knots must be a nondecreasing sequence of real numbers."
            );
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            var basisFunctions = new double[degree + 1];
            basisFunctions[0] = 1.0;

            var left = new double[degree + 1];
            var right = new double[degree + 1];

            for (int j = 1; j <= degree; j++)
            {
                left[j] = paramT - knots[spanIndex + 1 - j];
                right[j] = knots[spanIndex + j] - paramT;

                double saved = 0.0;

                for (int r = 0; r < j; r++)
                {
                    double temp = basisFunctions[r] / (right[r + 1] + left[j - r]);
                    basisFunctions[r] = saved + right[r + 1] * temp;
                    saved = left[j - r] * temp;
                }
                basisFunctions[j] = saved;
            }
            return basisFunctions;
        }

        /// <summary>
        /// Computes the non-vanishing basis functions and their derivatives up to the n-th derivative.
        /// </summary>
        /// <param name="spanIndex">The knot span index.</param>
        /// <param name="degree">The degree of the basis functions.</param>
        /// <param name="derivative">The order of derivatives to compute (k &lt;= degree).</param>
        /// <param name="knots">The knot vector.</param>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>Jagged array where result[k][j] is the k-th derivative of the j-th basis function.</returns>
        public static double[][] BasisFunctionsDerivatives(
            int spanIndex,
            int degree,
            int derivative,
            IReadOnlyList<double> knots,
            double paramT
        )
        {
            Validate.Argument(
                spanIndex >= 0,
                nameof(spanIndex),
                "SpanIndex must be greater than or equal zero."
            );
            Validate.Argument(
                degree >= 0,
                nameof(degree),
                "Degree must be greater than or equal zero."
            );
            Validate.Argument(
                derivative <= degree,
                nameof(derivative),
                "Derivative must not be greater than degree."
            );
            Validate.Argument(
                knots.Count > 0,
                nameof(knots),
                "Knots size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knots),
                nameof(knots),
                "Knots must be a nondecreasing sequence of real numbers."
            );
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            // Initialize jagged arrays
            var derivatives = new double[derivative + 1][];
            for (int i = 0; i <= derivative; i++)
                derivatives[i] = new double[degree + 1];

            var ndu = new double[degree + 1][];
            for (int i = 0; i <= degree; i++)
                ndu[i] = new double[degree + 1];

            ndu[0][0] = 1.0;

            var left = new double[degree + 1];
            var right = new double[degree + 1];

            double saved = 0.0;
            double temp = 0.0;

            for (int j = 1; j <= degree; j++)
            {
                left[j] = paramT - knots[spanIndex + 1 - j];
                right[j] = knots[spanIndex + j] - paramT;

                saved = 0.0;
                for (int r = 0; r < j; r++)
                {
                    ndu[j][r] = right[r + 1] + left[j - r];
                    temp = ndu[r][j - 1] / ndu[j][r];

                    ndu[r][j] = saved + right[r + 1] * temp;
                    saved = left[j - r] * temp;
                }
                ndu[j][j] = saved;
            }

            for (int j = 0; j <= degree; j++)
            {
                derivatives[0][j] = ndu[j][degree];
            }

            var a = new double[2][];
            a[0] = new double[degree + 1];
            a[1] = new double[degree + 1];

            for (int r = 0; r <= degree; r++)
            {
                int s1 = 0;
                int s2 = 1;
                a[0][0] = 1.0;

                for (int k = 1; k <= derivative; k++)
                {
                    double d = 0.0;
                    int rk = r - k;
                    int pk = degree - k;

                    if (r >= k)
                    {
                        a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                        d = a[s2][0] * ndu[rk][pk];
                    }

                    int j1 = 0;
                    int j2 = 0;

                    if (rk >= -1)
                        j1 = 1;
                    else
                        j1 = -rk;

                    if (r - 1 <= pk)
                        j2 = k - 1;
                    else
                        j2 = degree - r;

                    for (int j = j1; j <= j2; j++)
                    {
                        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                        d += a[s2][j] * ndu[rk + j][pk];
                    }
                    if (r <= pk)
                    {
                        a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                        d += a[s2][k] * ndu[r][pk];
                    }
                    derivatives[k][r] = d;

                    int swap = s1;
                    s1 = s2;
                    s2 = swap;
                }
            }

            int rVal = degree;
            for (int k = 1; k <= derivative; k++)
            {
                for (int j = 0; j <= degree; j++)
                {
                    derivatives[k][j] *= rVal;
                }
                rVal *= degree - k;
            }
            return derivatives;
        }

        /// <summary>
        /// Computes basis functions and their first derivatives (optimized).
        /// </summary>
        /// <param name="spanIndex">The knot span index.</param>
        /// <param name="degree">The degree.</param>
        /// <param name="knots">The knot vector.</param>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>Jagged array [2][p+1] containing basis functions and their first derivatives.</returns>
        public static double[][] BasisFunctionsFirstOrderDerivative(
            int spanIndex,
            int degree,
            IReadOnlyList<double> knots,
            double paramT
        )
        {
            Validate.Argument(
                spanIndex >= 0,
                nameof(spanIndex),
                "SpanIndex must be greater than or equal zero."
            );
            Validate.Argument(
                degree >= 0,
                nameof(degree),
                "Degree must be greater than or equal zero."
            );
            Validate.Argument(
                1 <= degree,
                nameof(degree),
                "Derivative must not be greater than degree."
            );
            Validate.Argument(
                knots.Count > 0,
                nameof(knots),
                "Knots size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knots),
                nameof(knots),
                "Knots must be a nondecreasing sequence of real numbers."
            );
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            var derivatives = new double[2][];
            derivatives[0] = new double[degree + 1];
            derivatives[1] = new double[degree + 1];

            var ndu = new double[degree + 1][];
            for (int i = 0; i <= degree; i++)
                ndu[i] = new double[degree + 1];

            ndu[0][0] = 1.0;

            var left = new double[degree + 1];
            var right = new double[degree + 1];

            double saved = 0.0;
            double temp = 0.0;

            for (int j = 1; j <= degree; j++)
            {
                left[j] = paramT - knots[spanIndex + 1 - j];
                right[j] = knots[spanIndex + j] - paramT;

                saved = 0.0;
                for (int r = 0; r < j; r++)
                {
                    ndu[j][r] = right[r + 1] + left[j - r];
                    temp = ndu[r][j - 1] / ndu[j][r];

                    ndu[r][j] = saved + right[r + 1] * temp;
                    saved = left[j - r] * temp;
                }
                ndu[j][j] = saved;
            }

            for (int j = 0; j <= degree; j++)
            {
                derivatives[0][j] = ndu[j][degree];
            }

            var a = new double[2][];
            a[0] = new double[degree + 1];
            a[1] = new double[degree + 1];

            for (int r = 0; r <= degree; r++)
            {
                int s1 = 0;
                int s2 = 1;
                a[0][0] = 1.0;

                const int k = 1;
                {
                    double d = 0.0;
                    int rk = r - k;
                    int pk = degree - k;

                    if (r >= k)
                    {
                        a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                        d = a[s2][0] * ndu[rk][pk];
                    }

                    int j1 = 0;
                    int j2 = 0;

                    if (rk >= -1)
                        j1 = 1;
                    else
                        j1 = -rk;

                    if (r - 1 <= pk)
                        j2 = k - 1;
                    else
                        j2 = degree - r;

                    for (int j = j1; j <= j2; j++)
                    {
                        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                        d += a[s2][j] * ndu[rk + j][pk];
                    }
                    if (r <= pk)
                    {
                        a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                        d += a[s2][k] * ndu[r][pk];
                    }
                    derivatives[k][r] = d;

                    int swap = s1;
                    s1 = s2;
                    s2 = swap;
                }
            }

            int rVal = degree;

            for (int j = 0; j <= degree; j++)
            {
                derivatives[1][j] *= rVal;
            }
            return derivatives;
        }

        /// <summary>
        /// Computes a single basis function N_{i,p}(t).
        /// </summary>
        /// <param name="spanIndex">The index i of the basis function.</param>
        /// <param name="degree">The degree p.</param>
        /// <param name="knots">The knot vector.</param>
        /// <param name="paramT">The parameter t.</param>
        /// <returns>The value of the basis function.</returns>
        public static double OneBasisFunction(
            int spanIndex,
            int degree,
            IReadOnlyList<double> knots,
            double paramT
        )
        {
            Validate.Argument(
                spanIndex >= 0,
                nameof(spanIndex),
                "SpanIndex must be greater than or equal zero."
            );
            Validate.Argument(degree > 0, nameof(degree), "Degree must be greater than zero.");
            Validate.Argument(
                knots.Count > 0,
                nameof(knots),
                "Knots size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knots),
                nameof(knots),
                "Knots must be a nondecreasing sequence of real numbers."
            );
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            int m = knots.Count - 1;
            if (
                (spanIndex == 0 && MathUtils.IsAlmostEqualTo(paramT, knots[0]))
                || (spanIndex == m - degree - 1 && MathUtils.IsAlmostEqualTo(paramT, knots[m]))
            )
            {
                return 1.0;
            }
            if (
                MathUtils.IsLessThan(paramT, knots[spanIndex])
                || MathUtils.IsGreaterThanOrEqual(paramT, knots[spanIndex + degree + 1])
            )
            {
                return 0.0;
            }

            var N = new double[degree + spanIndex + 1];
            for (int j = 0; j <= degree; j++)
            {
                if (
                    MathUtils.IsGreaterThanOrEqual(paramT, knots[spanIndex + j])
                    && MathUtils.IsLessThan(paramT, knots[spanIndex + j + 1])
                )
                {
                    N[j] = 1.0;
                }
            }

            for (int k = 1; k <= degree; k++)
            {
                double saved = 0.0;
                if (!MathUtils.IsAlmostEqualTo(N[0], 0.0))
                {
                    saved =
                        ((paramT - knots[spanIndex]) * N[0])
                        / (knots[spanIndex + k] - knots[spanIndex]);
                }

                for (int j = 0; j < degree - k + 1; j++)
                {
                    double knotLeft = knots[spanIndex + j + 1];
                    double knotRight = knots[spanIndex + j + k + 1];
                    if (MathUtils.IsAlmostEqualTo(N[j + 1], 0.0))
                    {
                        N[j] = saved;
                        saved = 0.0;
                    }
                    else
                    {
                        double temp = N[j + 1] / (knotRight - knotLeft);
                        N[j] = saved + (knotRight - paramT) * temp;
                        saved = (paramT - knotLeft) * temp;
                    }
                }
            }
            return N[0];
        }

        /// <summary>
        /// Computes the derivatives of a single basis function N_{i,p}(t).
        /// </summary>
        /// <param name="spanIndex">The index i.</param>
        /// <param name="degree">The degree p.</param>
        /// <param name="derivative">The number of derivatives to compute.</param>
        /// <param name="knots">The knot vector.</param>
        /// <param name="paramT">The parameter t.</param>
        /// <returns>Array containing the function value and its derivatives.</returns>
        public static double[] OneBasisFunctionDerivative(
            int spanIndex,
            int degree,
            int derivative,
            IReadOnlyList<double> knots,
            double paramT
        )
        {
            Validate.Argument(
                spanIndex >= 0,
                nameof(spanIndex),
                "SpanIndex must be greater than or equal zero."
            );
            Validate.Argument(degree > 0, nameof(degree), "Degree must be greater than zero.");
            Validate.Argument(
                derivative <= degree,
                nameof(derivative),
                "Derivative must not be greater than degree."
            );
            Validate.Argument(
                knots.Count > 0,
                nameof(knots),
                "Knots size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knots),
                nameof(knots),
                "Knots must be a nondecreasing sequence of real numbers."
            );
            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            var derivatives = new double[derivative + 1];

            if (
                MathUtils.IsLessThan(paramT, knots[spanIndex])
                || MathUtils.IsGreaterThanOrEqual(paramT, knots[spanIndex + degree + 1])
            )
            {
                return derivatives;
            }

            var N = new double[degree + 1][];
            for (int i = 0; i <= degree; i++)
                N[i] = new double[degree + 1];

            for (int j = 0; j <= degree; j++)
            {
                if (
                    MathUtils.IsGreaterThanOrEqual(paramT, knots[spanIndex + j])
                    && MathUtils.IsLessThan(paramT, knots[spanIndex + j + 1])
                )
                {
                    N[j][0] = 1.0;
                }
            }
            for (int k = 1; k <= degree; k++)
            {
                double saved = 0.0;
                if (!MathUtils.IsAlmostEqualTo(N[0][k - 1], 0.0))
                {
                    saved =
                        ((paramT - knots[spanIndex]) * N[0][k - 1])
                        / (knots[spanIndex + k] - knots[spanIndex]);
                }
                for (int j = 0; j < degree - k + 1; j++)
                {
                    double knotLeft = knots[spanIndex + j + 1];
                    double knotRight = knots[spanIndex + j + k + 1];

                    if (MathUtils.IsAlmostEqualTo(N[j + 1][k - 1], 0.0))
                    {
                        N[j][k] = saved;
                        saved = 0.0;
                    }
                    else
                    {
                        double temp = N[j + 1][k - 1] / (knotRight - knotLeft);
                        N[j][k] = saved + (knotRight - paramT) * temp;
                        saved = (paramT - knotLeft) * temp;
                    }
                }
            }

            derivatives[0] = N[0][degree];

            for (int k = 1; k <= derivative; k++)
            {
                var ND = new double[k + 1];
                for (int j = 0; j <= k; j++)
                {
                    ND[j] = N[j][degree - k];
                }
                for (int jj = 1; jj <= k; jj++)
                {
                    double saved = MathUtils.IsAlmostEqualTo(ND[0], 0.0)
                        ? 0.0
                        : ND[0] / (knots[spanIndex + degree - k + jj] - knots[spanIndex]);
                    for (int j = 0; j < k - jj + 1; j++)
                    {
                        double knotLeft = knots[spanIndex + j + 1];
                        double knotRight = knots[spanIndex + j + degree - k + jj + 1];

                        if (MathUtils.IsAlmostEqualTo(ND[j + 1], 0.0))
                        {
                            ND[j] = (degree - k + jj) * saved;
                            saved = 0.0;
                        }
                        else
                        {
                            double temp = ND[j + 1] / (knotRight - knotLeft);
                            ND[j] = (degree - k + jj) * (saved - temp);
                            saved = temp;
                        }
                    }
                }
                derivatives[k] = ND[0];
            }
            return derivatives;
        }

        /// <summary>
        /// Computes all non-zero basis functions of all degrees from 0 to p at the given span.
        /// </summary>
        /// <param name="spanIndex">The knot span index.</param>
        /// <param name="degree">The max degree p.</param>
        /// <param name="knots">The knot vector.</param>
        /// <param name="knot">The parameter/knot value.</param>
        /// <returns>Jagged array where result[j][r] is the value of basis function of degree j.</returns>
        public static double[][] AllBasisFunctions(
            int spanIndex,
            int degree,
            IReadOnlyList<double> knots,
            double knot
        )
        {
            Validate.Argument(
                spanIndex >= 0,
                nameof(spanIndex),
                "SpanIndex must be greater than or equal zero."
            );
            Validate.Argument(
                degree >= 0,
                nameof(degree),
                "Degree must be greater than or equal zero."
            );
            Validate.Argument(
                knots.Count > 0,
                nameof(knots),
                "Knots size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knots),
                nameof(knots),
                "Knots must be a nondecreasing sequence of real numbers."
            );
            Validate.Range(knot, knots[0], knots[knots.Count - 1], nameof(knot));

            var bases = new double[degree + 1][];
            for (int i = 0; i <= degree; i++)
                bases[i] = new double[degree + 1];

            bases[0][0] = 1.0;

            var left = new double[degree + 1];
            var right = new double[degree + 1];

            for (int j = 1; j <= degree; j++)
            {
                left[j] = knot - knots[spanIndex + 1 - j];
                right[j] = knots[spanIndex + j] - knot;
                double saved = 0.0;
                for (int r = 0; r < j; r++)
                {
                    bases[j][r] = right[r + 1] + left[j - r];
                    double temp = bases[r][j - 1] / bases[j][r];

                    bases[r][j] = saved + right[r + 1] * temp;
                    saved = left[j - r] * temp;
                }
                bases[j][j] = saved;
            }
            return bases;
        }

        /// <summary>
        /// Generates the change-of-basis matrix from Bezier to Power basis.
        /// </summary>
        /// <param name="degree">The degree.</param>
        /// <returns>The conversion matrix.</returns>
        public static double[][] BezierToPowerMatrix(int degree)
        {
            var matrix = new double[degree + 1][];
            for (int i = 0; i <= degree; i++)
                matrix[i] = new double[degree + 1];

            matrix[0][0] = matrix[degree][degree] = 1.0;
            if (degree % 2 != 0)
            {
                matrix[degree][0] = -1.0;
            }
            else
            {
                matrix[degree][0] = 1.0;
            }

            double sign = -1.0;
            for (int i = 1; i < degree; i++)
            {
                matrix[i][i] = MathUtils.Binomial(degree, i);
                matrix[i][0] = matrix[degree][degree - i] = sign * matrix[i][i];
                sign = -sign;
            }

            int k1 = (degree + 1) / 2;
            int pk = degree - 1;
            for (int k = 1; k < k1; k++)
            {
                sign = -1.0;
                for (int j = k + 1; j <= pk; j++)
                {
                    matrix[j][k] = matrix[pk][degree - j] =
                        sign
                        * MathUtils.Binomial(degree, k)
                        * MathUtils.Binomial(degree - k, j - k);
                    sign = -sign;
                }
                pk = pk - 1;
            }
            return matrix;
        }

        /// <summary>
        /// Generates the change-of-basis matrix from Power to Bezier basis.
        /// </summary>
        /// <param name="degree">The degree.</param>
        /// <param name="matrix">The Bezier-to-Power matrix to invert/convert.</param>
        /// <returns>The Power-to-Bezier conversion matrix.</returns>
        public static double[][] PowerToBezierMatrix(int degree, double[][] matrix)
        {
            var inverseMatrix = new double[degree + 1][];
            for (int i = 0; i <= degree; i++)
                inverseMatrix[i] = new double[degree + 1];

            for (int i = 0; i <= degree; i++)
            {
                inverseMatrix[i][0] = inverseMatrix[degree][i] = 1.0;
                inverseMatrix[i][i] = 1.0 / (matrix[i][i]);
            }

            int k1 = (degree + 1) / 2;
            int pk = degree - 1;

            for (int k = 1; k < k1; k++)
            {
                for (int j = k + 1; j <= pk; j++)
                {
                    double d = 0.0;
                    for (int i = k; i < j; i++)
                    {
                        d = d - matrix[j][i] * inverseMatrix[i][k];
                    }
                    inverseMatrix[j][k] = d / (matrix[j][j]);
                    inverseMatrix[pk][degree - j] = inverseMatrix[j][k];
                }
                pk = pk - 1;
            }
            return inverseMatrix;
        }
    }
}
