namespace Nurbsy
{
    public static class Validate
    {
        public static void Argument(bool isInvalid, string nameOf, string message)
        {
            if (isInvalid)
                throw new ArgumentOutOfRangeException(nameOf, message);
        }

        public static void Range<T>(T value, T min, T max, string nameOf)
            where T : IComparable<T>
        {
            if (value.CompareTo(min) < 0 || value.CompareTo(max) > 0)
                throw new ArgumentOutOfRangeException(
                    nameOf,
                    $"Value {value} is out of range [{min}, {max}]."
                );
        }

        public static bool IsValidBezier(int degree, int controlPointsCount)
        {
            return controlPointsCount == degree + 1;
        }

        public static bool IsValidKnots(IReadOnlyList<double> knots)
        {
            for (int i = 0; i < knots.Count - 1; i++)
            {
                // Verify non-decreasing order: knots[i] <= knots[i+1]
                if (knots[i] > knots[i + 1])
                    return false;
            }
            return true;
        }

        public static bool IsValidBSpline(int degree, int controlPointsCount, int knotsCount)
        {
            return (knotsCount - 1) == (controlPointsCount - 1) + degree + 1;
        }

        public static bool IsValidNURBS(int degree, int controlPointsCount, int knotsCount)
        {
            return (knotsCount - 1) == (controlPointsCount - 1) + degree + 1;
        }
    }
}
