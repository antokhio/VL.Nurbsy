using System.Collections.Immutable;
using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core;
using VL.Core.Import;

namespace VL.Nurbsy
{
    public abstract class BezierCurveNodeXD<T> : BezierCurveNode<T>
        where T : struct
    {
        protected static readonly int DefaultDegree = 1;

        private Optional<int> _degree;

        protected BezierCurveNodeXD(BezierCurve<T> curve)
            : base(curve) { }

        public abstract void SetControlPoints(IReadOnlyList<T> controlPoints);

        public void SetDegree(Optional<int> degree)
        {
            if (_degree != degree)
            {
                if (degree.HasValue)
                {
                    Degree = Math.Max(1, degree.Value);
                }
                else
                {
                    Degree = null;
                }
                _degree = degree;
            }
        }
    }

    [ProcessNode(Name = "BezierCurve (2D)")]
    public class BezierCurveNode2D : BezierCurveNodeXD<Vector2>
    {
        protected static readonly IReadOnlyList<Vector2> DefaultControlPoints =
        [
            new(-0.5f, 0f),
            new(0.5f, 0f),
        ];
        private IReadOnlyList<Vector2> _controlPoints;

        public BezierCurveNode2D()
            : base(new(DefaultDegree, DefaultControlPoints)) { }

        public override void SetControlPoints(IReadOnlyList<Vector2> controlPoints)
        {
            if (_controlPoints != controlPoints)
            {
                var cps = controlPoints?.Any() ?? false ? controlPoints : DefaultControlPoints;
                ControlPoints = cps.Select(cp => new ControlPoint<Vector2>(cp)).ToImmutableArray();

                _controlPoints = controlPoints;
            }
        }
    }

    [ProcessNode(Name = "BezierCurve (3D)")]
    public class BezierCurveNode3D : BezierCurveNodeXD<Vector3>
    {
        protected static readonly IReadOnlyList<Vector3> DefaultControlPoints =
        [
            new(-0.5f, 0f, 0f),
            new(0.5f, 0f, 0f),
        ];
        private IReadOnlyList<Vector3> _controlPoints;

        public BezierCurveNode3D()
            : base(new(DefaultDegree, DefaultControlPoints)) { }

        public override void SetControlPoints(IReadOnlyList<Vector3> controlPoints)
        {
            if (_controlPoints != controlPoints)
            {
                var cps = controlPoints?.Any() ?? false ? controlPoints : DefaultControlPoints;
                ControlPoints = cps.Select(cp => new ControlPoint<Vector3>(cp)).ToImmutableArray();
                _controlPoints = controlPoints;
            }
        }
    }
}
