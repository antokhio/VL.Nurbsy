using Stride.Core.Mathematics;

namespace Nurbsy.Algorithm
{
    public static class Intersection2D
    {
        public static CurveCurveIntersectionType ComputeRays(
            Vector2 point0,
            Vector2 vector0,
            Vector2 point1,
            Vector2 vector1,
            out double param0,
            out double param1,
            out Vector2 intersectPoint
        )
        {
            if (MathUtils.IsAlmostEqualTo(vector0, vector1))
            {
                if (MathUtils.IsAlmostEqualTo(point0, point1))
                {
                    intersectPoint = point0;
                    param0 = param1 = 0;
                    return CurveCurveIntersectionType.Intersecting;
                }
            }

            Validate.Argument(
                MathUtils.IsZero(vector0),
                nameof(vector0),
                "Vector0 must not be zero vector."
            );
            Validate.Argument(
                MathUtils.IsZero(vector1),
                nameof(vector1),
                "Vector1 must not be zero vector."
            );

            // 2D Cross Product (Determinant)
            double cross = vector0.X * vector1.Y - vector0.Y * vector1.X;
            Vector2 diff = point1 - point0;
            double coinCross = diff.X * vector1.Y - diff.Y * vector1.X;

            if (MathUtils.IsZero(cross))
            {
                if (MathUtils.IsZero(coinCross))
                {
                    intersectPoint = Vector2.Zero;
                    param0 = param1 = 0;
                    return CurveCurveIntersectionType.Coincident;
                }
                else
                {
                    intersectPoint = Vector2.Zero;
                    param0 = param1 = 0;
                    return CurveCurveIntersectionType.Parallel;
                }
            }

            // Cramer's Rule
            // t = (diff X v1) / (v0 X v1)
            double pd1Dot = diff.X * vector1.Y - diff.Y * vector1.X;
            param0 = pd1Dot / cross;

            // u = (diff X v0) / (v0 X v1)
            double pd2Dot = diff.X * vector0.Y - diff.Y * vector0.X;
            param1 = pd2Dot / cross;

            Vector2 rayP0 = point0 + vector0 * (float)param0;

            // In 2D, non-parallel lines always intersect.
            intersectPoint = rayP0;
            return CurveCurveIntersectionType.Intersecting;
        }
    }
}
