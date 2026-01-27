using Nurbsy;
using Stride.Core.Mathematics;
using Stride.Rendering.ProceduralModels;
using VL.Core.Import;
using VL.Nurbsy.Helpers;
using VL.Nurbsy.Rendering;

namespace VL.Nurbsy
{
    /// <summary>
    /// 2D Bezier surface mesh generator.
    /// Excludes normals computation.
    /// </summary>
    [ProcessNode(Name = "BezierSurfaceMesh (2D)")]
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
    [ProcessNode(Name = "BezierSurfaceMesh (3D)")]
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
