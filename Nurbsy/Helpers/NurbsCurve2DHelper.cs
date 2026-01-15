using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    internal static class NurbsCurve2DHelper
    {
        public static Vector2 GetPointOnCurve(NurbsCurve<Vector2> curve, double paramT)
        {
            var degree = curve.Degree;
            var knots = curve.Knots;
            var controlPoints = curve.ControlPoints;

            Validate.Range(paramT, knots[0], knots[knots.Count - 1], nameof(paramT));

            if (MathUtils.IsAlmostEqualTo(paramT, knots[0]))
            {
                return controlPoints[0].Value;
            }

            int n = controlPoints.Count - 1;
            if (MathUtils.IsAlmostEqualTo(paramT, knots[n + degree + 1]))
            {
                return controlPoints[n].Value;
            }

            int knotSpanIndex = Polynomials.GetKnotSpanIndex(degree, knots, paramT);
            int originMultiplicity = Polynomials.GetKnotMultiplicity(knots, paramT);

            int times = degree - originMultiplicity;

            // working in Homogeneous coordinates (w*x, w*y, w)
            var temp = new Vector3[times + 1];

            for (int i = 0; i <= times; i++)
            {
                var cp = controlPoints[knotSpanIndex - degree + i];
                double w = cp.Weight;
                temp[i] = new Vector3(cp.Value * (float)w, (float)w);
            }

            for (int j = 1; j <= times; j++)
            {
                for (int i = 0; i <= times - j; i++)
                {
                    double num = paramT - knots[knotSpanIndex - degree + j + i];
                    double den =
                        knots[i + knotSpanIndex + 1] - knots[knotSpanIndex - degree + j + i];
                    double alpha = num / den;

                    temp[i] = (float)alpha * temp[i + 1] + (float)(1.0 - alpha) * temp[i];
                }
            }

            Vector3 result = temp[0];
            if (MathUtils.IsZero(result.Z))
                return Vector2.Zero;

            return new Vector2(result.X, result.Y) / result.Z;
        }
    }
}
