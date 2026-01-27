using Stride.Core.Mathematics;
using Stride.Graphics;

namespace Nurbsy.Rendering.Helpers
{
    public static class NurbsCurveModelFactory
    {
        /// <summary>
        /// Generates a 3D tube mesh around a NURBS curve with flat caps.
        /// </summary>
        /// <param name="nurbsCurve">The 3D NURBS curve.</param>
        /// <param name="radius">The radius of the tube.</param>
        /// <param name="tesselation">The number of segments along the length of the curve.</param>
        /// <param name="segments">The number of radial segments (cross-section resolution).</param>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateCappedTube(
            NurbsCurve<Vector3> nurbsCurve,
            float radius,
            int tesselation,
            int segments
        )
        {
            var tes = Math.Max(2, tesselation);
            var seg = Math.Max(3, segments);

            // Body: (tesselation + 1) rings * (segments + 1) vertices (for UV seam)
            int numRings = tesselation + 1;
            int vertsPerRing = segments + 1;
            int bodyVertexCount = numRings * vertsPerRing;
            int bodyIndexCount = tesselation * segments * 6;

            // Caps: 1 center + (segments + 1) rim vertices
            int capVertexCount = 1 + (segments + 1);
            int capIndexCount = segments * 3;

            int totalVertexCount = bodyVertexCount + (2 * capVertexCount);
            int totalIndexCount = bodyIndexCount + (2 * capIndexCount);

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(
                totalVertexCount
            );
            var indices = GC.AllocateUninitializedArray<int>(totalIndexCount);

            int vOffset = 0;
            int iOffset = 0;

            double uStart = nurbsCurve.Knots[nurbsCurve.Degree];
            double uEnd = nurbsCurve.Knots[nurbsCurve.Knots.Count - 1 - nurbsCurve.Degree];
            double uRange = uEnd - uStart;

            // Store first and last frame for caps
            Vector3 startCenter = default,
                startTangent = default,
                startNormal = default,
                startBinormal = default;
            Vector3 endCenter = default,
                endTangent = default,
                endNormal = default,
                endBinormal = default;

            for (int i = 0; i < numRings; i++)
            {
                float t = i / (float)tesselation; // Normalized 0..1
                double u = uStart + t * uRange;

                // 1. Calculate Frame
                // We use the helper to get derivatives for frame calculation
                var derivatives = nurbsCurve.GetDerivatives(2, (float)u);
                Vector3 P = derivatives[0];
                Vector3 d1 = derivatives[1];
                Vector3 d2 = derivatives[2];

                Vector3 T = Vector3.Normalize(d1);
                // Frenet-Serret frame calculation
                Vector3 B = Vector3.Normalize(Vector3.Cross(T, d2));
                Vector3 N = Vector3.Normalize(Vector3.Cross(B, T));

                // Handle degenerate cases (straight lines where d2 is zero or parallel to d1)
                if (float.IsNaN(B.X) || B.LengthSquared() < 1e-6f)
                {
                    // Arbitrary axis method
                    Vector3 up =
                        Math.Abs(Vector3.Dot(T, Vector3.UnitY)) > 0.9f
                            ? Vector3.UnitZ
                            : Vector3.UnitY;
                    B = Vector3.Normalize(Vector3.Cross(T, up));
                    N = Vector3.Normalize(Vector3.Cross(B, T));
                }

                // Store frames for caps
                if (i == 0)
                {
                    startCenter = P;
                    startTangent = T;
                    startNormal = N;
                    startBinormal = B;
                }
                if (i == numRings - 1)
                {
                    endCenter = P;
                    endTangent = T;
                    endNormal = N;
                    endBinormal = B;
                }

                // 2. Generate Ring Vertices
                for (int j = 0; j <= segments; j++)
                {
                    float angle = (j / (float)segments) * MathUtil.TwoPi;
                    float cos = (float)Math.Cos(angle);
                    float sin = (float)Math.Sin(angle);

                    // Position: P + Radius * (cos*N + sin*B)
                    // Note: We use N and B to define the cross-section plane
                    Vector3 offset = radius * (cos * N + sin * B);
                    Vector3 position = P + offset;
                    Vector3 normal = Vector3.Normalize(offset); // Normal points radially out

                    vertices[vOffset + j] = new VertexPositionNormalTexture
                    {
                        Position = position,
                        Normal = normal,
                        TextureCoordinate = new Vector2((float)j / segments, t), // U wraps around, V goes along length
                    };
                }

                vOffset += vertsPerRing;
            }

            // 3. Generate Body Indices
            // Grid triangulation: [i, j], [i+1, j], [i, j+1] ...
            int ringOffset = 0;
            for (int i = 0; i < tesselation; i++)
            {
                int nextRingOffset = ringOffset + vertsPerRing;
                for (int j = 0; j < segments; j++)
                {
                    int current = ringOffset + j;
                    int next = current + 1;
                    int above = nextRingOffset + j;
                    int aboveNext = above + 1;

                    // Triangle 1
                    indices[iOffset++] = current;
                    indices[iOffset++] = above;
                    indices[iOffset++] = next;

                    // Triangle 2
                    indices[iOffset++] = next;
                    indices[iOffset++] = above;
                    indices[iOffset++] = aboveNext;
                }
                ringOffset += vertsPerRing;
            }

            GenerateFlatTubeCap(
                vertices,
                indices,
                ref vOffset,
                ref iOffset,
                startCenter,
                -startTangent,
                startNormal,
                startBinormal,
                radius,
                segments,
                isStart: true
            );

            GenerateFlatTubeCap(
                vertices,
                indices,
                ref vOffset,
                ref iOffset,
                endCenter,
                endTangent,
                endNormal,
                endBinormal,
                radius,
                segments,
                isStart: false
            );

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }

        /// <summary>
        /// Generates a flat cap geometry.
        /// </summary>
        private static void GenerateFlatTubeCap(
            VertexPositionNormalTexture[] vertices,
            int[] indices,
            ref int vOffset,
            ref int iOffset,
            Vector3 center,
            Vector3 normal,
            Vector3 axisU, // Usually Normal of the curve frame
            Vector3 axisV, // Usually Binormal of the curve frame
            float radius,
            int segments,
            bool isStart
        )
        {
            int centerIndex = vOffset;

            // 1. Center Vertex
            vertices[vOffset++] = new VertexPositionNormalTexture
            {
                Position = center,
                Normal = normal,
                TextureCoordinate = new Vector2(0.5f, 0.5f),
            };

            // 2. Rim Vertices
            // We duplicate these even though they match the tube positions because
            // the Normal needs to point in the direction of the tangent (flat cap), not radially.
            int rimStartIndex = vOffset;
            for (int j = 0; j <= segments; j++)
            {
                float angle = (j / (float)segments) * MathUtil.TwoPi;
                float cos = (float)Math.Cos(angle);
                float sin = (float)Math.Sin(angle);

                Vector3 position = center + radius * (cos * axisU + sin * axisV);

                // Simple planar UV mapping centered at 0.5, 0.5
                Vector2 uv = new Vector2(0.5f + 0.5f * cos, 0.5f + 0.5f * sin);

                vertices[vOffset++] = new VertexPositionNormalTexture
                {
                    Position = position,
                    Normal = normal,
                    TextureCoordinate = uv,
                };
            }

            // 3. Indices (Triangle CCW)
            for (int j = 0; j < segments; j++)
            {
                if (isStart)
                {
                    indices[iOffset++] = centerIndex;
                    indices[iOffset++] = rimStartIndex + j;
                    indices[iOffset++] = rimStartIndex + j + 1;
                }
                else
                {
                    indices[iOffset++] = centerIndex;
                    indices[iOffset++] = rimStartIndex + j + 1;
                    indices[iOffset++] = rimStartIndex + j;
                }
            }
        }
    }
}
