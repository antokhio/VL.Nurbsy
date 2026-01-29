using System.Collections.Immutable;
using Nurbsy;
using Nurbsy.Algorithm;
using Nurbsy.Helpers;
using Stride.Core.Mathematics;
using VL.Core.Import;
using VL.Model;

namespace VL.Nurbsy
{
    /// <summary>
    /// Additional base class for <see cref="NurbsSurfaceNode{T}"/>
    /// that ressolves Int2 Degree and Knots
    /// </summary>
    /// <typeparam name="T">Vector2, Vector3</typeparam>
    [ProcessNode]
    public abstract class NurbsSurfaceXDNode<T> : NurbsSurfaceNode<T>
        where T : struct
    {
        protected static readonly IReadOnlyList<double> DefaultKnotsU = [0.0, 0.0, 1.0, 1.0];
        protected static readonly IReadOnlyList<double> DefaultKnotsV = [0.0, 0.0, 1.0, 1.0];

        private Int2 _degree = default;

        protected NurbsSurfaceXDNode(NurbsSurface<T> surface)
            : base(surface) { }

        public abstract void SetControlPoints(IReadOnlyList<IReadOnlyList<T>> controlPoints);

        public void SetDegree(Int2 degree = default)
        {
            if (_degree != degree)
            {
                if (degree == default)
                {
                    DegreeU = DefaultDegreeU;
                    DegreeV = DefaultDegreeV;
                }
                else
                {
                    // TODO: Needs soft throw here
                    DegreeU = Math.Max(1, degree.X);
                    DegreeV = Math.Max(1, degree.Y);
                }

                _degree = degree;
            }
        }

        public void SetKnotsU(
            [Pin(Visibility = PinVisibility.Optional)] IReadOnlyList<double> knotsU
        )
        {
            if (KnotsU != knotsU)
            {
                KnotsU = knotsU;
            }
        }

        public void SetKnotsV(
            [Pin(Visibility = PinVisibility.Optional)] IReadOnlyList<double> knotsV
        )
        {
            if (KnotsV != knotsV)
            {
                KnotsV = knotsV;
            }
        }
    }

    [ProcessNode(Name = "NurbsSurface (2d)")]
    public class NurbsSurfaceNode2D : NurbsSurfaceXDNode<Vector2>
    {
        protected static readonly IReadOnlyList<IReadOnlyList<Vector2>> DefaultControlPoints =
        [
            [new(-0.5f, 0.5f), new(0.5f, 0.5f)],
            [new(-0.5f, -0.5f), new(0.5f, -0.5f)],
        ];

        private IReadOnlyList<IReadOnlyList<Vector2>> _controlPoints;

        public NurbsSurfaceNode2D()
            : base(
                new(
                    DefaultDegreeU,
                    DefaultDegreeV,
                    DefaultControlPoints,
                    DefaultKnotsU,
                    DefaultKnotsV
                )
            ) { }

        public override void SetControlPoints(IReadOnlyList<IReadOnlyList<Vector2>> controlPoints)
        {
            if (_controlPoints != controlPoints)
            {
                var cps = controlPoints?.Any() ?? false ? controlPoints : DefaultControlPoints;

                ControlPoints = cps.Select(row =>
                        (IReadOnlyList<ControlPoint<Vector2>>)
                            row.Select(cp => new ControlPoint<Vector2>(cp)).ToImmutableArray()
                    )
                    .ToImmutableArray();

                _controlPoints = controlPoints;
            }
        }

        public static NurbsSurface<Vector2> CreateDefault()
        {
            return new NurbsSurface<Vector2>(
                DefaultDegreeU,
                DefaultDegreeV,
                DefaultControlPoints,
                DefaultKnotsU,
                DefaultKnotsV
            );
        }
    }

    [ProcessNode(Name = "NurbsSurface (3d)")]
    public class NurbsSurfaceNode3D : NurbsSurfaceXDNode<Vector3>
    {
        protected static readonly IReadOnlyList<IReadOnlyList<Vector3>> DefaultControlPoints =
        [
            [new(-0.5f, 0.5f, 0f), new(0.5f, 0.5f, 0f)],
            [new(-0.5f, -0.5f, 0f), new(0.5f, -0.5f, 0f)],
        ];
        private IReadOnlyList<IReadOnlyList<Vector3>> _controlPoints;

        public NurbsSurfaceNode3D()
            : base(
                new(
                    DefaultDegreeU,
                    DefaultDegreeV,
                    DefaultControlPoints,
                    DefaultKnotsU,
                    DefaultKnotsV
                )
            ) { }

        public override void SetControlPoints(IReadOnlyList<IReadOnlyList<Vector3>> controlPoints)
        {
            if (_controlPoints != controlPoints)
            {
                var cps = controlPoints?.Any() ?? false ? controlPoints : DefaultControlPoints;

                ControlPoints = cps.Select(row =>
                        (IReadOnlyList<ControlPoint<Vector3>>)
                            row.Select(cp => new ControlPoint<Vector3>(cp)).ToImmutableArray()
                    )
                    .ToImmutableArray();

                _controlPoints = controlPoints;
            }
        }

        public static NurbsSurface<Vector3> CreateDefault()
        {
            return new NurbsSurface<Vector3>(
                DefaultDegreeU,
                DefaultDegreeV,
                DefaultControlPoints,
                DefaultKnotsU,
                DefaultKnotsV
            );
        }
    }

    [ProcessNode(Name = "NurbsSurface (3D Ruled)")]
    public class NurbsSurfaceRuledNode : NurbsSurfaceNode<Vector3>
    {
        protected static readonly NurbsCurve<Vector3> DefaultCurve0 = new(
            1,
            [new Vector3(-0.5f, 0.5f, 0f), new Vector3(-0.5f, 0.5f, 0f)],
            KnotsUtils.GenerateClampedKnots(1, 2)
        );
        protected static readonly NurbsCurve<Vector3> DefaultCurve1 = new(
            1,
            [new Vector3(-0.5f, -0.5f, 0f), new Vector3(-0.5f, -0.5f, 0f)],
            KnotsUtils.GenerateClampedKnots(1, 2)
        );

        private NurbsCurve<Vector3> _curve0;
        private NurbsCurve<Vector3> _curve1;

        public NurbsSurfaceRuledNode()
            : base(NurbsSurfaceHelper.CreateRulledSurface(DefaultCurve0, DefaultCurve1)) { }

        public void SetCurve0(NurbsCurve<Vector3> curve0)
        {
            if (_curve0 != curve0)
            {
                _curve0 = curve0;
                Invalidate();
            }
        }

        public void SetCurve1(NurbsCurve<Vector3> curve1)
        {
            if (_curve1 != curve1)
            {
                _curve1 = curve1;
                Invalidate();
            }
        }

        protected override void Build()
        {
            Output = NurbsSurfaceHelper.CreateRulledSurface(_curve0, _curve1);
        }
    }
}
