using Nurbsy.Algorithm;
using Stride.Core.Mathematics;
using Stride.Graphics;

namespace Nurbsy.Rendering.Helpers
{
    public static class NurbsSurfaceModelFactory
    {
        /// <summary>
        /// NURBS Surface 2D mesh generation, without normals calculation.
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshData(
            NurbsSurface<Vector2> surface,
            Int2 tesselation
        )
        {
            var tX = Math.Max(2, tesselation.X);
            var tY = Math.Max(2, tesselation.Y);

            // Determine parametric domain from knot vectors
            var uKnots = surface.KnotsU;
            var vKnots = surface.KnotsV;
            var uDegree = surface.DegreeU;
            var vDegree = surface.DegreeV;

            float uMin = (float)uKnots[uDegree];
            float uMax = (float)uKnots[uKnots.Count - 1 - uDegree];
            float uRange = uMax - uMin;

            float vMin = (float)vKnots[vDegree];
            float vMax = (float)vKnots[vKnots.Count - 1 - vDegree];
            float vRange = vMax - vMin;

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            for (int j = 0; j < tY; j++)
            {
                float vFactor = j / (float)(tY - 1);
                float v = vMin + vFactor * vRange;

                for (int i = 0; i < tX; i++)
                {
                    float uFactor = i / (float)(tX - 1);
                    float u = uMin + uFactor * uRange;

                    var point = surface.GetPointOnSurface(new Vector2(u, v));
                    vertices[j * tX + i] = new VertexPositionNormalTexture
                    {
                        Position = new Vector3(point, 0),
                        Normal = Vector3.UnitZ,
                        TextureCoordinate = new Vector2(uFactor, vFactor),
                    };
                }
            }

            // Stride uses CCW for Front Faces by default.
            var indices = Triangulation.GenerateGridIndicesCCW(tX, tY);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }

        /// <summary>
        /// NURBS Surface 3D mesh generation, with normals calculation.
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshData(
            NurbsSurface<Vector3> surface,
            Int2 tesselation
        )
        {
            var tX = Math.Max(2, tesselation.X);
            var tY = Math.Max(2, tesselation.Y);

            // Determine parametric domain from knot vectors
            var uKnots = surface.KnotsU;
            var vKnots = surface.KnotsV;
            var uDegree = surface.DegreeU;
            var vDegree = surface.DegreeV;

            // Valid domain: [knots[degree], knots[count - 1 - degree]]
            float uMin = (float)uKnots[uDegree];
            float uMax = (float)uKnots[uKnots.Count - 1 - uDegree];
            float uRange = uMax - uMin;

            float vMin = (float)vKnots[vDegree];
            float vMax = (float)vKnots[vKnots.Count - 1 - vDegree];
            float vRange = vMax - vMin;

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            for (int j = 0; j < tY; j++)
            {
                float vFactor = j / (float)(tY - 1);
                float v = vMin + vFactor * vRange;

                for (int i = 0; i < tX; i++)
                {
                    float uFactor = i / (float)(tX - 1);
                    float u = uMin + uFactor * uRange;

                    // Get point and first derivatives efficiently
                    surface.GetRationalFirstOrderDerivatives(
                        new Vector2(u, v),
                        out var point,
                        out var su,
                        out var sv
                    );

                    // Compute normal via cross product of partial derivatives
                    // Reversed (sv x su) to match CW winding
                    var normal = Vector3.Cross(su, sv);
                    if (normal.LengthSquared() > 1e-8f)
                    {
                        normal = Vector3.Normalize(normal);
                    }
                    else
                    {
                        normal = Vector3.UnitZ; // Degenerate handling
                    }

                    vertices[j * tX + i] = new VertexPositionNormalTexture
                    {
                        Position = point,
                        Normal = normal,
                        TextureCoordinate = new Vector2(uFactor, vFactor),
                    };
                }
            }

            // 3D Surfaces use CW indices to correct culling
            var indices = Triangulation.GenerateGridIndicesCW(tX, tY);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }
    }
}
