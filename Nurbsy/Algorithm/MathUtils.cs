namespace Nurbsy.Algorithm
{
    internal static class MathUtils
    {
        public const double Epsilon = 1e-9;

        public static bool IsAlmostEqualTo(double a, double b, double epsilon = Epsilon)
        {
            return Math.Abs(a - b) <= epsilon;
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
