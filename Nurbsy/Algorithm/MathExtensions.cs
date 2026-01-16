using Stride.Core.Mathematics;

namespace Nurbsy.Algorithm
{
    public static class MathExtensions
    {
        /// <summary>
        /// Creates a random vector orthogonal to the specified vector.
        /// </summary>
        /// <param name="xyz">The vector to be orthogonal to.</param>
        /// <returns>A normalized random orthogonal vector.</returns>
        public static Vector3 RandomOrthogonal(this ref Vector3 xyz)
        {
            // 1. Normalize the input vector
            Vector3 normal = Vector3.Normalize(xyz);

            // 2. Create the "swizzled" vector based on the C++ snippet: (-normal.Z, normal.X, normal.Y)
            Vector3 swizzled = new Vector3(-normal.Z, normal.X, normal.Y);

            // 3. Calculate tangent: normal x swizzled
            Vector3.Cross(ref normal, ref swizzled, out Vector3 tangent);

            // 4. Calculate bitangent: normal x tangent
            Vector3.Cross(ref normal, ref tangent, out Vector3 bitangent);

            // 5. Generate random angle between -Pi and Pi
            // Using System.Random.Shared (thread-safe, .NET 6+)
            float angle = System.Random.Shared.NextSingle() * MathF.PI * 2.0f - MathF.PI;

            // 6. Calculate the result: tangent * sin(angle) + bitangent * cos(angle)
            Vector3 result = (tangent * MathF.Sin(angle)) + (bitangent * MathF.Cos(angle));

            // 7. Normalize and return
            result.Normalize();
            return result;
        }

        /// <summary>
        /// Creates a random vector orthogonal to the specified vector.
        /// </summary>
        /// <param name="vector">The vector to be orthogonal to.</param>
        /// <returns>A normalized orthogonal vector (randomly rotated 90 or -90 degrees).</returns>
        public static Vector2 RandomOrthogonal(this ref Vector2 vector)
        {
            // 1. Normalize the input to ensure we return a unit vector
            Vector2 normal = Vector2.Normalize(vector);

            // 2. Create the perpendicular vector (-y, x)
            // This is mathematically equivalent to a 90-degree rotation
            Vector2 orthogonal = new Vector2(-normal.Y, normal.X);

            // 3. Randomly flip the direction to pick one of the two possible orthogonal vectors
            // (Either 90 degrees or -90 degrees relative to the original)
            if (System.Random.Shared.Next(2) == 0)
            {
                orthogonal.X = -orthogonal.X;
                orthogonal.Y = -orthogonal.Y;
            }

            return orthogonal;
        }
    }
}
