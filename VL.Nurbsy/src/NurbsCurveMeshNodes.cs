using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
    [ProcessNode(Name = "NurbsCurveMesh (3D)")]
    public class NurbsCurveMeshNode3D : NurbsCurveMeshNode<Vector3>
    {
        public NurbsCurveMeshNode3D()
            : base(NurbsCurveNode3D.CreateDefault()) { }

        public override void SetCurve(NurbsCurve<Vector3> curve)
        {
            base.SetCurve(curve);
        }
    }
}
