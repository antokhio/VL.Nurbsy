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
            var (tX, tY) = GetTessellation(tesselation);
            var cps = surface.ControlPoints;
            var uParams = TextureParametrization.Resolve(surface.DegreeU, surface.KnotsU, cps.Count, tX);
            var vParams = TextureParametrization.Resolve(surface.DegreeV, surface.KnotsV, cps[0].Count, tY);
            return Generate2D(surface, tX, tY, uParams, vParams);
        }

        /// <summary>
        /// NURBS Surface 3D mesh generation, with normals calculation.
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshData(
            NurbsSurface<Vector3> surface,
            Int2 tesselation
        )
        {
            var (tX, tY) = GetTessellation(tesselation);
            var cps = surface.ControlPoints;
            var uParams = TextureParametrization.Resolve(surface.DegreeU, surface.KnotsU, cps.Count, tX);
            var vParams = TextureParametrization.Resolve(surface.DegreeV, surface.KnotsV, cps[0].Count, tY);
            return Generate3D(surface, tX, tY, uParams, vParams);
        }

        /// <summary>
        /// Cached 2D mesh generation. Identical output to
        /// <see cref="GenerateMeshData(NurbsSurface{Vector2}, Int2)"/>, but reuses the
        /// tessellation parameters from <paramref name="cache"/> while only control point
        /// positions change (see <see cref="NurbsSurfaceParamCache{T}"/>).
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshDataCached(
            NurbsSurface<Vector2> surface,
            Int2 tesselation,
            NurbsSurfaceParamCache<Vector2> cache
        )
        {
            var (tX, tY) = GetTessellation(tesselation);
            cache.GetOrBuild(in surface, tX, tY);
            return Generate2D(surface, tX, tY, cache.UParams, cache.VParams);
        }

        /// <summary>
        /// Cached 3D mesh generation. Identical output to
        /// <see cref="GenerateMeshData(NurbsSurface{Vector3}, Int2)"/>, but reuses the
        /// tessellation parameters from <paramref name="cache"/> while only control point
        /// positions change (see <see cref="NurbsSurfaceParamCache{T}"/>).
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshDataCached(
            NurbsSurface<Vector3> surface,
            Int2 tesselation,
            NurbsSurfaceParamCache<Vector3> cache
        )
        {
            var (tX, tY) = GetTessellation(tesselation);
            cache.GetOrBuild(in surface, tX, tY);
            return Generate3D(surface, tX, tY, cache.UParams, cache.VParams);
        }

        private static (int X, int Y) GetTessellation(Int2 tesselation)
        {
            return (Math.Max(2, tesselation.X), Math.Max(2, tesselation.Y));
        }

        // Texture coordinates are the evenly spaced tessellation factors, while the surface
        // parameters come from TextureParametrization (topology-only). The texture is therefore
        // bound to the control net: moving a control point moves its part of the texture.
        private static GeometricMeshData<VertexPositionNormalTexture> Generate2D(
            NurbsSurface<Vector2> surface,
            int tX,
            int tY,
            float[] uParams,
            float[] vParams
        )
        {
            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            for (int i = 0; i < tX; i++)
            {
                float uFactor = i / (float)(tX - 1);
                float u = uParams[i];

                for (int j = 0; j < tY; j++)
                {
                    float vFactor = j / (float)(tY - 1);
                    float v = vParams[j];

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

        private static GeometricMeshData<VertexPositionNormalTexture> Generate3D(
            NurbsSurface<Vector3> surface,
            int tX,
            int tY,
            float[] uParams,
            float[] vParams
        )
        {
            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            for (int i = 0; i < tX; i++)
            {
                float uFactor = i / (float)(tX - 1);
                float u = uParams[i];

                for (int j = 0; j < tY; j++)
                {
                    float vFactor = j / (float)(tY - 1);
                    float v = vParams[j];

                    surface.GetRationalFirstOrderDerivatives(
                        new Vector2(u, v),
                        out var point,
                        out var su,
                        out var sv
                    );

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