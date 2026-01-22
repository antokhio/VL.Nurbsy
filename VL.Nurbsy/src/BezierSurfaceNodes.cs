using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
    [ProcessNode(Name = "BezierSurface (2D)")]
    public class BezierSurface2DNode : BezierSurfaceNode<Vector2>
    {
        static readonly IReadOnlyList<IReadOnlyList<Vector2>> DefaultControlPoints =
        [
            [new Vector2(-0.5f, 0.5f), new Vector2(0.5f, 0.5f)],
            [new Vector2(-0.5f, -0.5f), new Vector2(0.5f, -0.5f)],
        ];

        public BezierSurface2DNode()
            : base(new(DefaultControlPointsCount, DefaultControlPoints)) { }

        public override void SetControlPoints(IReadOnlyList<IReadOnlyList<Vector2>> controlPoints)
        {
            base.SetControlPoints(controlPoints);
        }
    }
}
