using System.Collections.Immutable;
using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
    [ProcessNode(Name = "BezierSurface (2D)")]
    public class BezierSurfaceNode2D : BezierSurfaceNode<Vector2>
    {
        internal static readonly IReadOnlyList<IReadOnlyList<Vector2>> DefaultControlPoints =
        [
            [new(-0.5f, 0.5f), new(0.5f, 0.5f)],
            [new(-0.5f, -0.5f), new(0.5f, -0.5f)],
        ];

        private IReadOnlyList<IReadOnlyList<Vector2>> _controlPoints;

        public BezierSurfaceNode2D()
            : base(new(DefaultDegreeU, DefaultDegreeV, DefaultControlPoints)) { }

        public void SetControlPoints(IReadOnlyList<IReadOnlyList<Vector2>> controlPoints)
        {
            if (controlPoints != _controlPoints)
            {
                var cps = controlPoints?.Any() ?? false ? controlPoints : DefaultControlPoints;

                DegreeU = cps.Count - 1;
                DegreeV = cps[0].Count - 1;

                ControlPoints = cps.Select(row =>
                        (IReadOnlyList<ControlPoint<Vector2>>)
                            row.Select(cp => new ControlPoint<Vector2>(cp)).ToImmutableArray()
                    )
                    .ToImmutableArray();

                _controlPoints = controlPoints;
            }
        }
    }

    [ProcessNode(Name = "BezierSurface (3D)")]
    public class BezierSurfaceNode3D : BezierSurfaceNode<Vector3>
    {
        internal static readonly IReadOnlyList<IReadOnlyList<Vector3>> DefaultControlPoints =
        [
            [new(-0.5f, 0.5f, 0f), new(0.5f, 0.5f, 0f)],
            [new(-0.5f, -0.5f, 0f), new(0.5f, -0.5f, 0f)],
        ];

        private IReadOnlyList<IReadOnlyList<Vector3>> _controlPoints;

        public BezierSurfaceNode3D()
            : base(new(DefaultDegreeU, DefaultDegreeV, DefaultControlPoints)) { }

        public void SetControlPoints(IReadOnlyList<IReadOnlyList<Vector3>> controlPoints)
        {
            if (controlPoints != _controlPoints)
            {
                var cps = controlPoints?.Any() ?? false ? controlPoints : DefaultControlPoints;

                var degreeU = cps.Count + 1;
                var degreeV = cps[0].Count + 1;

                ControlPoints = cps.Select(row =>
                        (IReadOnlyList<ControlPoint<Vector3>>)
                            row.Select(cp => new ControlPoint<Vector3>(cp)).ToImmutableArray()
                    )
                    .ToImmutableArray();

                _controlPoints = controlPoints;
            }
        }
    }
}
