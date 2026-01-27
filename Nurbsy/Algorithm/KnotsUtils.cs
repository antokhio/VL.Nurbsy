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
    public static class KnotsUtils
    {
        public static IReadOnlyList<double> GenerateClampedKnots(int degree, int controlPointCount)
        {
            int knotCount = controlPointCount + degree + 1;
            List<double> knots = new List<double>(knotCount);
            // Add clamped start knots
            for (int i = 0; i <= degree; i++)
            {
                knots.Add(0.0);
            }
            // Add internal knots
            int internalKnotCount = knotCount - 2 * (degree + 1);
            for (int i = 1; i <= internalKnotCount; i++)
            {
                knots.Add((double)i / (internalKnotCount + 1));
            }
            // Add clamped end knots
            for (int i = 0; i <= degree; i++)
            {
                knots.Add(1.0);
            }
            return knots;
        }

        public static int GetContinuity(int degree, IReadOnlyList<double> knotVector, double knot)
        {
            Validate.Argument(degree > 0, nameof(degree), "Degree must be greater than zero.");
            int multi = Polynomials.GetKnotMultiplicity(knotVector, knot);
            return degree - multi;
        }

        public static IReadOnlyList<double> Rescale(
            IReadOnlyList<double> knotVector,
            double min,
            double max
        )
        {
            double originMin = knotVector[0];
            double originMax = knotVector[knotVector.Count - 1];
            double k = (max - min) / (originMax - originMin);

            int size = knotVector.Count;
            var result = new double[size];
            for (int i = 0; i < size; i++)
            {
                result[i] = k * (knotVector[i] - originMin) + min;
            }
            return result;
        }

        public static double GetNode(int degree, IReadOnlyList<double> knotVector, int index)
        {
            double sum = 0.0;
            for (int j = 1; j <= degree; j++)
            {
                sum += knotVector[index + j];
            }
            return sum / degree;
        }

        internal static Dictionary<double, int> GetKnotMultiplicityMap(
            IReadOnlyList<double> knotVector
        )
        {
            var map = new Dictionary<double, int>();
            foreach (var t in knotVector)
            {
                if (map.ContainsKey(t))
                    map[t]++;
                else
                    map[t] = 1;
            }
            return map;
        }

        internal static Dictionary<double, int> GetInternalKnotMultiplicityMap(
            IReadOnlyList<double> knotVector
        )
        {
            var map = new Dictionary<double, int>();
            // Usually internal knots are those strictly between first and last, but robust definition handles multiplicities
            // Assuming clamped knot vector, first degree+1 are min, last degree+1 are max.
            if (knotVector.Count == 0)
                return map;

            double first = knotVector[0];
            double last = knotVector[knotVector.Count - 1];

            foreach (var t in knotVector)
            {
                if (MathUtils.IsAlmostEqualTo(t, first) || MathUtils.IsAlmostEqualTo(t, last))
                    continue;

                if (map.ContainsKey(t))
                    map[t]++;
                else
                    map[t] = 1;
            }
            return map;
        }

        public static IReadOnlyList<double> AverageKnotVector(
            int degree,
            IReadOnlyList<double> parameters
        )
        {
            int n = parameters.Count;
            int m = n + degree + 1;
            var knots = new double[m];

            for (int i = 0; i <= degree; i++)
            {
                knots[i] = 0.0;
                knots[m - 1 - i] = 1.0;
            }

            for (int j = 1; j < n - degree; j++)
            {
                double sum = 0.0;
                for (int k = j; k < j + degree; k++)
                {
                    sum += parameters[k];
                }
                knots[j + degree] = sum / degree;
            }

            return knots;
        }

        public static bool IsClamped(int degree, IReadOnlyList<double> knotVector)
        {
            double first = knotVector[0];
            for (int i = 1; i <= degree; i++)
            {
                if (!MathUtils.IsAlmostEqualTo(knotVector[i], first))
                    return false;
            }

            double last = knotVector[knotVector.Count - 1];
            for (int i = knotVector.Count - 2; i >= knotVector.Count - 1 - degree; i--)
            {
                if (!MathUtils.IsAlmostEqualTo(knotVector[i], last))
                    return false;
            }

            return true;
        }

        public static bool IsUniform(IReadOnlyList<double> knotVector)
        {
            if (knotVector.Count < 2)
                return true;

            double delta = knotVector[1] - knotVector[0];
            for (int i = 1; i < knotVector.Count - 1; i++)
            {
                double d = knotVector[i + 1] - knotVector[i];
                if (!MathUtils.IsAlmostEqualTo(d, delta))
                    return false;
            }
            return true;
        }

        /// <summary>
        /// Reverse a knot vector while preserving spacing.
        /// </summary>
        /// <param name="knots">The original knot vector.</param>
        /// <returns>The reversed knot vector.</returns>
        public static IReadOnlyList<double> ReverseKnots(IReadOnlyList<double> knots)
        {
            int size = knots.Count;
            var reversed = new double[size];
            reversed[0] = knots[0];

            for (int i = 1; i < size; i++)
            {
                reversed[i] = reversed[i - 1] + (knots[size - i] - knots[size - i - 1]);
            }

            return reversed;
        }

        public static List<double> GetUniqueKnots(IReadOnlyList<double> knots)
        {
            var unique = new List<double> { knots[0] };
            for (int i = 1; i < knots.Count; i++)
            {
                if (!MathUtils.IsAlmostEqualTo(knots[i], unique[^1]))
                {
                    unique.Add(knots[i]);
                }
            }
            return unique;
        }

        public static void GetInsertionKnots(
            IReadOnlyList<double> knots1,
            IReadOnlyList<double> knots2,
            out List<double> insert1,
            out List<double> insert2
        )
        {
            insert1 = new List<double>();
            insert2 = new List<double>();

            // Collect all unique knot values from both vectors
            var allKnots = knots1.Concat(knots2).OrderBy(x => x).ToList();
            var uniqueKnots = new List<double>();

            if (allKnots.Count > 0)
            {
                uniqueKnots.Add(allKnots[0]);
                for (int i = 1; i < allKnots.Count; i++)
                {
                    if (!MathUtils.IsAlmostEqualTo(allKnots[i], uniqueKnots[uniqueKnots.Count - 1]))
                    {
                        uniqueKnots.Add(allKnots[i]);
                    }
                }
            }

            // Determine how many times we need to insert each unique knot to unify multiplicity
            foreach (var u in uniqueKnots)
            {
                int count1 = CountKnots(knots1, u);
                int count2 = CountKnots(knots2, u);
                int max = Math.Max(count1, count2);

                for (int i = 0; i < max - count1; i++)
                    insert1.Add(u);
                for (int i = 0; i < max - count2; i++)
                    insert2.Add(u);
            }
        }

        private static int CountKnots(IReadOnlyList<double> knots, double val)
        {
            int count = 0;
            foreach (var k in knots)
            {
                if (MathUtils.IsAlmostEqualTo(k, val))
                    count++;
            }
            return count;
        }
    }
}
