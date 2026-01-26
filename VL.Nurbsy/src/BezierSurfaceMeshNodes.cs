using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// 2D Bezier surface mesh generator.
    /// Excludes normals computation.
    /// </summary>
    [ProcessNode(Name = "BezierSurfaceMesh (2D)")]
    public class BezierSurfaceMeshNode2D : BezierSurfaceMeshNode<Vector2>
    {
        public BezierSurfaceMeshNode2D()
            : base(BezierSurfaceNode2D.CreateDefault()) { }

        public override void SetSurface(BezierSurface<Vector2> surface)
        {
            base.SetSurface(surface);
        }
    }

    /// <summary>
    /// 3D Bezier surface mesh generator.
    /// Includes normals computation.
    /// </summary>
    [ProcessNode(Name = "BezierSurfaceMesh (3D)")]
    public class BezierSurfaceMeshNode3D : BezierSurfaceMeshNode<Vector3>
    {
        public BezierSurfaceMeshNode3D()
            : base(BezierSurfaceNode3D.CreateDefault()) { }

        public override void SetSurface(BezierSurface<Vector3> surface)
        {
            base.SetSurface(surface);
        }
    }
}
