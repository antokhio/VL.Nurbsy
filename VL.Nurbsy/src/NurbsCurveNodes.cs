using System.Collections.Immutable;
using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// Base class for managing NurbsCurve instance, that definines
    /// uniformly weighted controlpoints as Vector2 or Vector3
    /// </summary>
    /// <typeparam name="T"></typeparam>
    [ProcessNode]
    public abstract class NurbsCurveNodeXD<T> : NurbsCurveNode<T>
        where T : struct
    {
        internal static readonly IReadOnlyList<double> DefaultKnots = [0.0, 0.0, 1.0, 1.0];

        protected NurbsCurveNodeXD(NurbsCurve<T> curve)
            : base(curve) { }

        public abstract void SetControlPoints(IReadOnlyList<T> controlPoints);

        public void SetDegree(int degree = DefaultDegree)
        {
            if (Degree != degree)
            {
                // TODO: soft throw here
                Degree = Math.Max(1, degree);
            }
        }

        public void SetKnots(
            [Pin(Visibility = Model.PinVisibility.Optional)] IReadOnlyList<double> knots
        )
        {
            if (knots != Knots)
            {
                Knots = knots;
            }
        }
    }

    [ProcessNode(Name = "NurbsCurve (2D)")]
    public class NurbsCurveNode2D : NurbsCurveNodeXD<Vector2>
    {
        static readonly IReadOnlyList<Vector2> DefaultControlPoints =
        [
            new(-0.5f, 0f),
            new(0.5f, 0f),
        ];

        private IReadOnlyList<Vector2> _controlPoints;

        public NurbsCurveNode2D()
            : base(new(DefaultDegree, DefaultControlPoints, DefaultKnots)) { }

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

    [ProcessNode(Name = "NurbsCurve (3D)")]
    public class NurbsCurveNode3D : NurbsCurveNodeXD<Vector3>
    {
        static readonly IReadOnlyList<Vector3> DefaultControlPoints =
        [
            new(-0.5f, 0f, 0f),
            new(0.5f, 0f, 0f),
        ];

        private IReadOnlyList<Vector3> _controlPoints;

        public NurbsCurveNode3D()
            : base(new(DefaultDegree, DefaultControlPoints, DefaultKnots)) { }

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
