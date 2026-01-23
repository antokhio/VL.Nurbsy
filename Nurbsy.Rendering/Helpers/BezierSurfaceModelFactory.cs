using Nurbsy.Algorithm;
using Stride.Core.Mathematics;
using Stride.Graphics;

namespace Nurbsy.Rendering.Helpers
{
    public static class BezierSurfaceModelFactory
    {
        public static GeometricMeshData<VertexPositionNormalTexture> Create(
            BezierSurface<Vector2> surface,
            Int2 tesselation
        )
        {
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
                    vertices[y * tX + x] = new VertexPositionNormalTexture
                    {
                        Position = new Vector3(point, 0),
                        Normal = Vector3.UnitZ,
                        TextureCoordinate = new Vector2(u, v),
                    };
                }
            }

            var indices = Triangulation.GenerateCCWGrid(tX, tY);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }
    }
}
