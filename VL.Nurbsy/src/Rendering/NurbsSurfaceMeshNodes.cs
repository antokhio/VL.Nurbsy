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
    /// Uses the exact arc-length parametrization (see <see cref="NurbsSurface{T}.Sample(Vector2)"/>)
    /// to avoid tessellation compression, recomputed in full whenever the surface changes. Heavier
    /// on CPU than the "(Cached) (Experimental)" variant, which is preferable when control points
    /// are animated every frame and this cost becomes a bottleneck.
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
    /// Uses the exact arc-length parametrization (see <see cref="NurbsSurface{T}.Sample(Vector2)"/>)
    /// to avoid tessellation compression, recomputed in full whenever the surface changes. Heavier
    /// on CPU than the "(Cached) (Experimental)" variant, which is preferable when control points
    /// are animated every frame and this cost becomes a bottleneck.
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
    /// Experimental, cached fast-path 2D Nurbs surface mesh generator.
    /// Excludes normals computation.
    /// Caches the arc-length lookup tables used to correct tessellation compression and only
    /// rebuilds them when the surface's degree or control point grid dimensions change - not on
    /// every control point position update. This makes it much cheaper for animated control
    /// points, at the cost of assuming the V-direction arc length is roughly constant across all
    /// U columns (see <see cref="NurbsSurfaceArcLengthCache{T}"/>).
    /// </summary>
    [ProcessNode(Name = "NurbsSurface (2d Mesh)")]
    public class NurbsSurfaceMeshNode2DCached : SurfaceMeshNode<NurbsSurface<Vector2>>
    {
        private readonly NurbsSurfaceArcLengthCache<Vector2> _cache = new();

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
    /// Experimental, cached fast-path 3D Nurbs surface mesh generator.
    /// Includes normals computation.
    /// Caches the arc-length lookup tables used to correct tessellation compression and only
    /// rebuilds them when the surface's degree or control point grid dimensions change - not on
    /// every control point position update. This makes it much cheaper for animated control
    /// points, at the cost of assuming the V-direction arc length is roughly constant across all
    /// U columns (see <see cref="NurbsSurfaceArcLengthCache{T}"/>).
    /// </summary>
    [ProcessNode(Name = "NurbsSurface (3d Mesh)")]
    public class NurbsSurfaceMeshNode3DCached : SurfaceMeshNode<NurbsSurface<Vector3>>
    {
        private readonly NurbsSurfaceArcLengthCache<Vector3> _cache = new();

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
