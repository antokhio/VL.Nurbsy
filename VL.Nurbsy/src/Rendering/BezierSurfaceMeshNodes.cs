using Nurbsy;
using Stride.Core.Mathematics;
using Stride.Rendering.ProceduralModels;
using VL.Core.Import;
using VL.Nurbsy.Helpers;

namespace VL.Nurbsy.Rendering
{
    /// <summary>
    /// 2D Bezier surface mesh generator.
    /// Excludes normals computation.
    /// </summary>
    [ProcessNode(Name = "BezierSurface (2d Mesh)")]
    public class BezierSurfaceMeshNode2D : SurfaceMeshNode<BezierSurface<Vector2>>
    {
        public BezierSurfaceMeshNode2D()
            : base(BezierSurfaceNode2D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new BezierSurfaceModel<Vector2>(Surface, Tesselation);
        }

        protected override object GetResourceKey()
        {
            return (GetGameProvider(), typeof(BezierSurfaceModel<Vector2>), Surface, Tesselation);
        }
    }

    /// <summary>
    /// 3D Bezier surface mesh generator.
    /// Includes normals computation.
    /// </summary>
    [ProcessNode(Name = "BezierSurface (3d Mesh)")]
    public class BezierSurfaceMeshNode3D : SurfaceMeshNode<BezierSurface<Vector3>>
    {
        public BezierSurfaceMeshNode3D()
            : base(BezierSurfaceNode3D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new BezierSurfaceModel<Vector3>(Surface, Tesselation);
        }

        protected override object GetResourceKey()
        {
            return (GetGameProvider(), typeof(BezierSurfaceModel<Vector3>), Surface, Tesselation);
        }
    }
}
