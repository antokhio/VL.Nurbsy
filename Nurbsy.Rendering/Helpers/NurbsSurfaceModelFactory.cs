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

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            // Reference iso-curve along U, taken at the middle of the V domain, used to map
            // normalized tessellation factor tu -> real surface parameter u. Uniform parameter
            // spacing does not yield uniform geometric spacing for degree > 1, so we mirror the
            // arc-length correction used by NurbsSurface<T>.Sample to avoid tessellation compression.
            var knotsV = surface.KnotsV;
            float vRef = (float)((knotsV[0] + knotsV[knotsV.Count - 1]) * 0.5);
            var refCurveU = surface.GetCurveV(vRef);

            for (int i = 0; i < tX; i++)
            {
                float uFactor = i / (float)(tX - 1);
                float u = (float)refCurveU.GetParamAt(uFactor);

                // Iso-curve along V at the resolved u, used to map tv -> real surface parameter v.
                var curveV = surface.GetCurveU(u);

                for (int j = 0; j < tY; j++)
                {
                    float vFactor = j / (float)(tY - 1);
                    float v = (float)curveV.GetParamAt(vFactor);

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

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            // Reference iso-curve along U, taken at the middle of the V domain, used to map
            // normalized tessellation factor tu -> real surface parameter u. Uniform parameter
            // spacing does not yield uniform geometric spacing for degree > 1, so we mirror the
            // arc-length correction used by NurbsSurface<T>.Sample to avoid tessellation compression.
            var knotsV = surface.KnotsV;
            float vRef = (float)((knotsV[0] + knotsV[knotsV.Count - 1]) * 0.5);
            var refCurveU = surface.GetCurveV(vRef);

            for (int i = 0; i < tX; i++)
            {
                float uFactor = i / (float)(tX - 1);
                float u = (float)refCurveU.GetParamAt(uFactor);

                // Iso-curve along V at the resolved u, used to map tv -> real surface parameter v.
                var curveV = surface.GetCurveU(u);

                for (int j = 0; j < tY; j++)
                {
                    float vFactor = j / (float)(tY - 1);
                    float v = (float)curveV.GetParamAt(vFactor);

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

        /// <summary>
        /// Experimental, cached fast-path for 2D mesh generation. Uses precomputed arc-length
        /// lookup tables (see <see cref="NurbsSurfaceArcLengthCache{T}"/>) instead of resolving the
        /// exact arc-length parametrization for every vertex, at the cost of assuming the V-direction
        /// arc length is roughly constant across all U columns. Intended for animated control points,
        /// where the cache avoids re-running the (relatively expensive) arc-length integration every frame.
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshDataCached(
            NurbsSurface<Vector2> surface,
            Int2 tesselation,
            NurbsSurfaceArcLengthCache<Vector2> cache
        )
        {
            var tX = Math.Max(2, tesselation.X);
            var tY = Math.Max(2, tesselation.Y);

            cache.GetOrBuild(in surface);
            var uLut = cache.ULut;
            var vLut = cache.VLut;

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            for (int i = 0; i < tX; i++)
            {
                float uFactor = i / (float)(tX - 1);
                float u = (float)uLut.Sample(uFactor);

                for (int j = 0; j < tY; j++)
                {
                    float vFactor = j / (float)(tY - 1);
                    float v = (float)vLut.Sample(vFactor);

                    var point = surface.GetPointOnSurface(new Vector2(u, v));
                    vertices[j * tX + i] = new VertexPositionNormalTexture
                    {
                        Position = new Vector3(point, 0),
                        Normal = Vector3.UnitZ,
                        TextureCoordinate = new Vector2(uFactor, vFactor),
                    };
                }
            }

            var indices = Triangulation.GenerateGridIndicesCCW(tX, tY);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }

        /// <summary>
        /// Experimental, cached fast-path for 3D mesh generation. See the 2D overload's remarks;
        /// the same arc-length LUT caching and approximation applies here, with normals still
        /// computed exactly per-vertex via the resolved (u, v).
        /// </summary>
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateMeshDataCached(
            NurbsSurface<Vector3> surface,
            Int2 tesselation,
            NurbsSurfaceArcLengthCache<Vector3> cache
        )
        {
            var tX = Math.Max(2, tesselation.X);
            var tY = Math.Max(2, tesselation.Y);

            cache.GetOrBuild(in surface);
            var uLut = cache.ULut;
            var vLut = cache.VLut;

            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(tX * tY);

            for (int i = 0; i < tX; i++)
            {
                float uFactor = i / (float)(tX - 1);
                float u = (float)uLut.Sample(uFactor);

                for (int j = 0; j < tY; j++)
                {
                    float vFactor = j / (float)(tY - 1);
                    float v = (float)vLut.Sample(vFactor);

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

            var indices = Triangulation.GenerateGridIndicesCW(tX, tY);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }
    }
}
