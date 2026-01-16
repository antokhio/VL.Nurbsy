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
    }
}
