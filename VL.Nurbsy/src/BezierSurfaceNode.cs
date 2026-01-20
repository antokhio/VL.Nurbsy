using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
    [ProcessNode()]
    public abstract class BezierSurfaceNode<T>
        where T : struct
    {
        internal static readonly Int2 DefaultControlPointsCount = new(2, 2);

        private IReadOnlyList<IReadOnlyList<T>> _controlPoints;
        private Int2 _controlPointsCount = DefaultControlPointsCount;
        private bool _invalidate = false;

        public BezierSurface<T> Output { get; protected set; }

        protected BezierSurfaceNode(BezierSurface<T> surface)
        {
            Output = surface;
        }

        public virtual void SetControlPoints(IReadOnlyList<IReadOnlyList<T>> controlPoints)
        {
            if (_controlPoints != controlPoints)
            {
                _controlPoints = controlPoints;

                _invalidate = true;
            }
        }

        public virtual void SetControlPointsCount(Int2 controlPointsCount)
        {
            if (_controlPointsCount != controlPointsCount)
            {
                _controlPointsCount = controlPointsCount;

                _invalidate = true;
            }
        }

        public virtual void Update()
        {
            if (_invalidate)
            {
                Output = new BezierSurface<T>(_controlPointsCount, _controlPoints);

                _invalidate = false;
            }
        }
    }

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
