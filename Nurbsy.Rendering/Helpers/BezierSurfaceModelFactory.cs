using Nurbsy.Algorithm;
using Stride.Core.Mathematics;
using Stride.Graphics;

namespace Nurbsy.Rendering.Helpers
{
    public static class BezierSurfaceModelFactory
    {
        /// <summary>
        /// 2D mesh generation for flat surfaces, without normals computation.
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshData(
            BezierSurface<Vector2> surface,
            Int2 tesselation
        )
        {
            var tX = Math.Max(2, tesselation.X);
            var tY = Math.Max(2, tesselation.Y);

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            for (int j = 0; j < tY; j++)
            {
                for (int i = 0; i < tX; i++)
                {
                    var u = i / (float)(tX - 1);
                    var v = j / (float)(tY - 1);
                    var point = surface.GetPointOnSurface(new Vector2(u, v));
                    vertices[j * tX + i] = new VertexPositionNormalTexture
                    {
                        Position = new Vector3(point, 0),
                        Normal = Vector3.UnitZ,
                        TextureCoordinate = new Vector2(u, v),
                    };
                }
            }

            var indices = Triangulation.GenerateGridIndeciesCW(tX, tY);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }

        /// <summary>
        /// 3D mesh generation for flat surfaces, with normal computation.
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshData(
            BezierSurface<Vector3> surface,
            Int2 tesselation
        )
        {
            float eps = Constants.DistanceEpsilon;

            var tX = Math.Max(2, tesselation.X);
            var tY = Math.Max(2, tesselation.Y);

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            for (int y = 0; y < tY; y++)
            {
                for (int x = 0; x < tX; x++)
                {
                    var u = x / (float)(tX - 1);
                    var v = y / (float)(tY - 1);

                    var point = surface.GetPointOnSurface(new Vector2(u, v));

                    // Compute normals using finite differences

                    var u0 = Math.Max(0, u - eps);
                    var u1 = Math.Min(1, u + eps);

                    var v0 = Math.Max(0, v - eps);
                    var v1 = Math.Min(1, v + eps);

                    // dP/du ≈ (P(u+e) - P(u-e)) / 2e
                    var pLeft = surface.GetPointOnSurface(new Vector2(u0, v));
                    var pRight = surface.GetPointOnSurface(new Vector2(u1, v));
                    var tangentU = pRight - pLeft;

                    // dP/dv ≈ (P(v+e) - P(v-e)) / 2e
                    var pDown = surface.GetPointOnSurface(new Vector2(u, v0));
                    var pUp = surface.GetPointOnSurface(new Vector2(u, v1));
                    var tangentV = pUp - pDown;

                    var normal = Vector3.Cross(tangentU, tangentV);

                    if (normal.LengthSquared() > float.Epsilon)
                    {
                        normal.Normalize();
                    }
                    else
                    {
                        normal = Vector3.UnitZ;
                    }

                    vertices[y * tX + x] = new VertexPositionNormalTexture
                    {
                        Position = point,
                        Normal = normal,
                        TextureCoordinate = new Vector2(u, v),
                    };
                }
            }

            var indices = Triangulation.GenerateGridIndeciesCW(tX, tY);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }
    }
}
