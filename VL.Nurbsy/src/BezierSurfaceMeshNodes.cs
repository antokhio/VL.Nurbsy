using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
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
}
