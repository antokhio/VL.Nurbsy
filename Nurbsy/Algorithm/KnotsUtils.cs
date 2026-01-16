using System.Runtime.InteropServices;

namespace Nurbsy.Algorithm
{
    internal static class KnotsUtils
    {
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

        public static IReadOnlyList<double> GetInsertedKnotElement(
            int degree,
            IReadOnlyList<double> knotVector,
            double startParam,
            double endParam
        )
        {
            Validate.Argument(
                degree >= 0,
                nameof(degree),
                "Degree must be greater than or equal to zero."
            );
            Validate.Argument(
                knotVector.Count > 0,
                nameof(knotVector),
                "KnotVector size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(knotVector),
                nameof(knotVector),
                "KnotVector must be a nondecreasing sequence of real numbers."
            );

            var result = new List<double>();
            int startMulti = Polynomials.GetKnotMultiplicity(knotVector, startParam);
            if (startMulti < degree)
            {
                for (int i = 0; i < degree - startMulti; i++)
                {
                    result.Add(startParam);
                }
            }

            int endMulti = Polynomials.GetKnotMultiplicity(knotVector, endParam);
            if (endMulti < degree)
            {
                for (int i = 0; i < degree - endMulti; i++)
                {
                    result.Add(endParam);
                }
            }

            return result;
        }

        public static IReadOnlyDictionary<double, int> GetKnotMultiplicityMap(
            IReadOnlyList<double> knotVector
        )
        {
            var result = new Dictionary<double, int>();

            int count = knotVector.Count;
            for (int i = 0; i < count; i++)
            {
                double knot = knotVector[i];
                if (!result.ContainsKey(knot))
                {
                    int multi = Polynomials.GetKnotMultiplicity(knotVector, knot);
                    result.Add(knot, multi);
                }
            }
            return result;
        }

        public static IReadOnlyDictionary<double, int> GetInternalKnotMultiplicityMap(
            IReadOnlyList<double> knotVector
        )
        {
            // We need a mutable dictionary initially to modify it
            var result = new Dictionary<double, int>();

            // Re-implementing essentially what GetKnotMultiplicityMap does but mutable locally
            int count = knotVector.Count;
            for (int i = 0; i < count; i++)
            {
                double knot = knotVector[i];
                if (!result.ContainsKey(knot))
                {
                    int multi = Polynomials.GetKnotMultiplicity(knotVector, knot);
                    result.Add(knot, multi);
                }
            }

            if (result.Count > 0)
            {
                // Remove first (min)
                var minKey = result.Keys.Min();
                result.Remove(minKey);

                // Remove last (max) if it still exists
                if (result.Count > 0)
                {
                    var maxKey = result.Keys.Max();
                    result.Remove(maxKey);
                }
            }
            return result;
        }

        public static void GetInsertedKnotElement(
            IReadOnlyList<double> knotVector0,
            IReadOnlyList<double> knotVector1,
            List<double> insertElements0,
            List<double> insertElements1
        )
        {
            var map0 = GetKnotMultiplicityMap(knotVector0);
            var map1 = GetKnotMultiplicityMap(knotVector1);

            foreach (var kvp0 in map0)
            {
                double key0 = kvp0.Key;
                int count0 = kvp0.Value;

                if (!map1.TryGetValue(key0, out int count1))
                {
                    for (int i = 0; i < count0; i++)
                    {
                        insertElements1.Add(key0);
                    }
                }
                else
                {
                    if (count0 > count1)
                    {
                        int times = count0 - count1;
                        for (int j = 0; j < times; j++)
                        {
                            insertElements1.Add(key0);
                        }
                    }
                    else
                    {
                        int times = count1 - count0;
                        for (int j = 0; j < times; j++)
                        {
                            insertElements0.Add(key0);
                        }
                    }
                }
            }

            foreach (var kvp1 in map1)
            {
                double key1 = kvp1.Key;
                int count1 = kvp1.Value;

                if (!map0.ContainsKey(key1))
                {
                    for (int i = 0; i < count1; i++)
                    {
                        insertElements0.Add(key1);
                    }
                }
            }

            insertElements0.Sort();
            insertElements1.Sort();
        }

        public static IReadOnlyList<IReadOnlyList<double>> GetInsertedKnotElements(
            IReadOnlyList<IReadOnlyList<double>> knotVectors
        )
        {
            var finalMap = new Dictionary<double, int>();
            for (int i = 0; i < knotVectors.Count; i++)
            {
                var kv = knotVectors[i];
                var map = GetKnotMultiplicityMap(kv);
                foreach (var kvp in map)
                {
                    double key = kvp.Key;
                    int count = kvp.Value;
                    ref int currentVal = ref CollectionsMarshal.GetValueRefOrAddDefault(
                        finalMap,
                        key,
                        out bool exists
                    );
                    if (!exists)
                    {
                        currentVal = count;
                    }
                    else
                    {
                        if (currentVal < count)
                        {
                            currentVal = count;
                        }
                    }
                }
            }

            var result = new List<IReadOnlyList<double>>(knotVectors.Count);
            for (int i = 0; i < knotVectors.Count; i++)
            {
                var kv = knotVectors[i];
                var map = GetKnotMultiplicityMap(kv);

                var insertElements = new List<double>();
                foreach (var kvp in finalMap)
                {
                    double key = kvp.Key;
                    int count = kvp.Value;

                    if (!map.TryGetValue(key, out int currentCount))
                    {
                        for (int j = 0; j < count; j++)
                        {
                            insertElements.Add(key);
                        }
                    }
                    else
                    {
                        int times = count - currentCount;
                        for (int j = 0; j < times; j++)
                        {
                            insertElements.Add(key);
                        }
                    }
                }
                insertElements.Sort();
                result.Add(insertElements);
            }
            return result;
        }

        public static IReadOnlyList<double> GetInsertedKnotElements(
            int insertKnotsNumber,
            IReadOnlyList<double> knotVector
        )
        {
            var uniqueKnotVector = knotVector.Distinct().ToList();
            uniqueKnotVector.Sort();

            var insert = new List<double>(insertKnotsNumber);
            InsertMidKnotCore(uniqueKnotVector, insert, insertKnotsNumber);
            insert.Sort();
            return insert;
        }

        public static bool IsUniform(IReadOnlyList<double> knotVector)
        {
            var map = GetKnotMultiplicityMap(knotVector);
            if (map.Count == 0)
            {
                return false;
            }

            var knots = map.Keys.ToList();
            knots.Sort();

            if (map[knots[0]] != map[knots[knots.Count - 1]])
            {
                return false;
            }

            if (knots.Count < 2)
                return true;

            double standard = knots[1] - knots[0];
            for (int i = 1; i < knots.Count - 1; i++)
            {
                double current = knots[i];
                double next = knots[i + 1];

                double gap = next - current;
                if (!MathUtils.IsAlmostEqualTo(gap, standard))
                {
                    return false;
                }

                int currentCounts = map[current];
                int nextCounts = map[next];
                if (currentCounts != nextCounts)
                {
                    if (i + 1 == knots.Count - 1)
                    {
                        continue;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        private static void InsertMidKnotCore(
            List<double> uniqueKnotVector,
            List<double> insert,
            int limitNumber
        )
        {
            if (insert.Count == limitNumber)
            {
                return;
            }
            else
            {
                double standard = Constants.DoubleEpsilon;
                int index = -1;
                for (int i = 0; i < uniqueKnotVector.Count - 1; i++)
                {
                    double delta = uniqueKnotVector[i + 1] - uniqueKnotVector[i];
                    if (delta > standard)
                    {
                        standard = delta;
                        index = i;
                    }
                }

                if (index != -1)
                {
                    double current = uniqueKnotVector[index] + standard * 0.5;
                    uniqueKnotVector.Add(current);
                    uniqueKnotVector.Sort();
                    insert.Add(current);

                    InsertMidKnotCore(uniqueKnotVector, insert, limitNumber);
                }
            }
        }
    }
}
