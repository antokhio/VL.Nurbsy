namespace Nurbsy
{
    public static class NurbsCurveHelper
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
    }
}
