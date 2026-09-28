using Nurbsy;
using Nurbsy.Rendering;
using Nurbsy.Rendering.Helpers;
using Stride.Core.Mathematics;
using Stride.Rendering.ProceduralModels;
using VL.Core.Import;
using VL.Nurbsy.Helpers;

namespace VL.Nurbsy.Rendering
{
    /// <summary>
    /// 2D Nurbs surface mesh generator.
    /// Excludes normals computation.
    /// Texture coordinates follow the control net (see <see cref="TextureParametrization"/>).
    /// Recomputes the tessellation parameters on every change; prefer "NurbsSurface (2d Mesh)".
    /// </summary>
    [ProcessNode(Name = "NurbsSurface (2d Mesh Obsolete)")]
    public class NurbsSurfaceMeshNode2D : SurfaceMeshNode<NurbsSurface<Vector2>>
    {
        public NurbsSurfaceMeshNode2D()
            : base(NurbsSurfaceNode2D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsSurfaceModel<Vector2>(Surface, Tesselation);
        }

        protected override object GetResourceKey()
        {
            return (GetGameProvider(), typeof(NurbsSurfaceModel<Vector2>), Surface, Tesselation);
        }
    }

    /// <summary>
    /// 3D Nurbs surface mesh generator.
    /// Includes normals computation.
    /// Texture coordinates follow the control net (see <see cref="TextureParametrization"/>).
    /// Recomputes the tessellation parameters on every change; prefer "NurbsSurface (3d Mesh)".
    /// </summary>
    [ProcessNode(Name = "NurbsSurface (3d Mesh Obsolete)")]
    public class NurbsSurfaceMeshNode3D : SurfaceMeshNode<NurbsSurface<Vector3>>
    {
        public NurbsSurfaceMeshNode3D()
            : base(NurbsSurfaceNode3D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsSurfaceModel<Vector3>(Surface, Tesselation);
        }

        protected override object GetResourceKey()
        {
            return (GetGameProvider(), typeof(NurbsSurfaceModel<Vector3>), Surface, Tesselation);
        }
    }

    /// <summary>
    /// 2D Nurbs surface mesh generator.
    /// Excludes normals computation.
    /// Texture coordinates follow the control net (see <see cref="TextureParametrization"/>):
    /// moving a control point drags its part of the texture, while an evenly spaced grid gives an
    /// uncompressed texture. Tessellation parameters are cached and only rebuilt when degree,
    /// knots or grid dimensions change (see <see cref="NurbsSurfaceParamCache{T}"/>).
    /// </summary>
    [ProcessNode(Name = "NurbsSurface (2d Mesh)")]
    public class NurbsSurfaceMeshNode2DCached : SurfaceMeshNode<NurbsSurface<Vector2>>
    {
        private readonly NurbsSurfaceParamCache<Vector2> _cache = new();

        public NurbsSurfaceMeshNode2DCached()
            : base(NurbsSurfaceNode2D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsSurfaceModelCached<Vector2>(Surface, Tesselation, _cache);
        }

        protected override object GetResourceKey()
        {
            return (
                GetGameProvider(),
                typeof(NurbsSurfaceModelCached<Vector2>),
                Surface,
                Tesselation
            );
        }
    }

    /// <summary>
    /// 3D Nurbs surface mesh generator.
    /// Includes normals computation.
    /// Texture coordinates follow the control net (see <see cref="TextureParametrization"/>):
    /// moving a control point drags its part of the texture, while an evenly spaced grid gives an
    /// uncompressed texture. Tessellation parameters are cached and only rebuilt when degree,
    /// knots or grid dimensions change (see <see cref="NurbsSurfaceParamCache{T}"/>).
    /// </summary>
    [ProcessNode(Name = "NurbsSurface (3d Mesh)")]
    public class NurbsSurfaceMeshNode3DCached : SurfaceMeshNode<NurbsSurface<Vector3>>
    {
        private readonly NurbsSurfaceParamCache<Vector3> _cache = new();

        public NurbsSurfaceMeshNode3DCached()
            : base(NurbsSurfaceNode3D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsSurfaceModelCached<Vector3>(Surface, Tesselation, _cache);
        }

        protected override object GetResourceKey()
        {
            return (
                GetGameProvider(),
                typeof(NurbsSurfaceModelCached<Vector3>),
                Surface,
                Tesselation
            );
        }
    }
}
