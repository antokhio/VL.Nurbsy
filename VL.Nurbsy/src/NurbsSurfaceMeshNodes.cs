using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// 2D Bezier surface mesh generator.
    /// Excludes normals computation.
    /// </summary>
    [ProcessNode(Name = "NurbsSurfaceMesh (2D)")]
    public class NurbsSurfaceMeshNode2D : NurbsSurfaceMeshNode<Vector2>
    {
        public NurbsSurfaceMeshNode2D()
            : base(NurbsSurfaceNode2D.CreateDefault()) { }

        public override void SetSurface(NurbsSurface<Vector2> surface)
        {
            base.SetSurface(surface);
        }
    }

    /// <summary>
    /// 3D Bezier surface mesh generator.
    /// Includes normals computation.
    /// </summary>
    [ProcessNode(Name = "NurbsSurfaceMesh (3D)")]
    public class NurbsSurfaceMeshNode3D : NurbsSurfaceMeshNode<Vector3>
    {
        public NurbsSurfaceMeshNode3D()
            : base(NurbsSurfaceNode3D.CreateDefault()) { }

        public override void SetSurface(NurbsSurface<Vector3> surface)
        {
            base.SetSurface(surface);
        }
    }
}
