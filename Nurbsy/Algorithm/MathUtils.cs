using Stride.Core.Mathematics;

namespace Nurbsy.Algorithm
{
    internal static class MathUtils
    {
        public const double Epsilon = 1e-9;

        public static bool IsAlmostEqualTo(double a, double b, double epsilon = Epsilon)
        {
            return Math.Abs(a - b) <= epsilon;
        }

        public static bool IsAlmostEqualTo(float a, float b, double epsilon = Epsilon)
        {
            return Math.Abs(a - b) <= epsilon;
        }

        public static bool IsAlmostEqualTo(Vector2 a, Vector2 b, double epsilon = Epsilon)
        {
            return Vector2.DistanceSquared(a, b) <= (float)(epsilon * epsilon);
        }

        public static bool IsAlmostEqualTo(Vector3 a, Vector3 b, double epsilon = Epsilon)
        {
            return Vector3.DistanceSquared(a, b) <= (float)(epsilon * epsilon);
        }

        public static bool IsGreaterThanOrEqual(double a, double b, double epsilon = Epsilon)
        {
            return a >= b - epsilon;
        }

        public static bool IsLessThanOrEqual(double a, double b, double epsilon = Epsilon)
        {
            return a <= b + epsilon;
        }

        public static bool IsLessThan(double a, double b, double epsilon = Epsilon)
        {
            return a < b - epsilon;
        }

        public static bool IsGreaterThan(double a, double b, double epsilon = Epsilon)
        {
            return a > b + epsilon;
        }

        public static bool IsZero(double a, double epsilon = Epsilon)
        {
            return Math.Abs(a) <= epsilon;
        }

        public static bool IsZero(Vector2 v, double epsilon = Epsilon)
        {
            return v.LengthSquared() <= (float)(epsilon * epsilon);
        }

        public static bool IsZero(Vector3 v, double epsilon = Epsilon)
        {
            return v.LengthSquared() <= (float)(epsilon * epsilon);
        }

        public static double Binomial(int n, int k)
        {
            if (k < 0 || k > n)
                return 0.0;
            if (k == 0 || k == n)
                return 1.0;

            if (k > n / 2)
                k = n - k;

            double res = 1.0;
            for (int i = 1; i <= k; ++i)
            {
                res = res * (n - i + 1) / i;
            }
            return res;
        }
    }
}
