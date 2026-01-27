using Nurbsy;
using Nurbsy.Rendering;
using Stride.Core.Mathematics;
using Stride.Rendering.ProceduralModels;
using VL.Core.Import;
using VL.Nurbsy.Helpers;

namespace VL.Nurbsy
{
    [ProcessNode(Name = "NurbsCurveMesh (3D Tube)")]
    public class NurbsCurveTubeMeshNode3D : CurveTubeMeshNode<NurbsCurve<Vector3>>
    {
        public NurbsCurveTubeMeshNode3D()
            : base(NurbsCurveNode3D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsCurveModel<Vector3>(Curve, Radius, Tesselation, Segments);
        }

        protected override object GetResourceKey()
        {
            return (
                GetGameProvider(),
                typeof(NurbsCurveModel<Vector3>),
                Curve,
                Radius,
                Tesselation,
                Segments
            );
        }
    }
}
