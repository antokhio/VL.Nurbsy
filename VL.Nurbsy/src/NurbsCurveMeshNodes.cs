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
            return new NurbsCurveTubeModel<Vector3>(Curve, Radius, Tesselation, Segments);
        }

        protected override object GetResourceKey()
        {
            return (
                GetGameProvider(),
                typeof(NurbsCurveTubeModel<Vector3>),
                Curve,
                Radius,
                Tesselation,
                Segments
            );
        }
    }

    [ProcessNode(Name = "NurbsCurveMesh (2D Tape)")]
    public class NurbsCurveTapeMeshNode2D : CurveTapeMeshNode<NurbsCurve<Vector2>>
    {
        public NurbsCurveTapeMeshNode2D()
            : base(NurbsCurveNode2D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsCurveTapeModel<Vector2>(Curve, Width, Tesselation, Direction);
        }

        protected override object GetResourceKey()
        {
            return (
                GetGameProvider(),
                typeof(NurbsCurveTapeModel<Vector2>),
                Curve,
                Width,
                Tesselation,
                Direction
            );
        }
    }

    [ProcessNode(Name = "NurbsCurveMesh (3D Tape)")]
    public class NurbsCurveTapeMeshNode3D : CurveTapeMeshNode<NurbsCurve<Vector3>>
    {
        public NurbsCurveTapeMeshNode3D()
            : base(NurbsCurveNode3D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsCurveTapeModel<Vector3>(Curve, Width, Tesselation, Direction);
        }

        protected override object GetResourceKey()
        {
            return (
                GetGameProvider(),
                typeof(NurbsCurveTapeModel<Vector3>),
                Curve,
                Width,
                Tesselation,
                Direction
            );
        }
    }
}
