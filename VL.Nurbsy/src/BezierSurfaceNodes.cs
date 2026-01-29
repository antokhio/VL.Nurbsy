using System.Collections.Immutable;
using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// Additional base class for <see cref="BezierSurfaceNode{T}"/>
    /// that ressolves Int2 Degree
    /// </summary>
    /// <typeparam name="T">Vector2, Vector3</typeparam>
    [ProcessNode]
    public abstract class BezierSurfaceNodeXD<T> : BezierSurfaceNode<T>
        where T : struct
    {
        protected static readonly Int2 DefaultDegree = new(1, 1);

        private Optional<Int2> _degree;

        protected BezierSurfaceNodeXD(BezierSurface<T> surface)
            : base(surface) { }

        public abstract void SetControlPoints(IReadOnlyList<IReadOnlyList<T>> controlPoints);

        public void SetDegree(Optional<Int2> degree)
        {
            if (_degree != degree)
            {
                if (degree.HasValue)
                {
                    DegreeU = Math.Max(1, degree.Value.X);
                    DegreeV = Math.Max(1, degree.Value.Y);
                }
                else
                {
                    DegreeU = null;
                    DegreeV = null;
                }

                _degree = degree;
            }
        }
    }

    [ProcessNode(Name = "BezierSurface (2d)")]
    public class BezierSurfaceNode2D : BezierSurfaceNodeXD<Vector2>
    {
        protected static readonly IReadOnlyList<IReadOnlyList<Vector2>> DefaultControlPoints =
        [
            [new(-0.5f, 0.5f), new(0.5f, 0.5f)],
            [new(-0.5f, -0.5f), new(0.5f, -0.5f)],
        ];

        private IReadOnlyList<IReadOnlyList<Vector2>> _controlPoints;

        public BezierSurfaceNode2D()
            : base(BezierSurfaceNode2D.CreateDefault()) { }

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

        public static BezierSurface<Vector2> CreateDefault() =>
            new(DefaultDegree.X, DefaultDegree.Y, DefaultControlPoints);
    }

    [ProcessNode(Name = "BezierSurface (3d)")]
    public class BezierSurfaceNode3D : BezierSurfaceNodeXD<Vector3>
    {
        protected static readonly IReadOnlyList<IReadOnlyList<Vector3>> DefaultControlPoints =
        [
            [new(-0.5f, 0.5f, 0f), new(0.5f, 0.5f, 0f)],
            [new(-0.5f, -0.5f, 0f), new(0.5f, -0.5f, 0f)],
        ];

        private IReadOnlyList<IReadOnlyList<Vector3>> _controlPoints;

        public BezierSurfaceNode3D()
            : base(new(DefaultDegree.X, DefaultDegree.Y, DefaultControlPoints)) { }

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

        public static BezierSurface<Vector3> CreateDefault() =>
            new(DefaultDegree.X, DefaultDegree.Y, DefaultControlPoints);
    }
}
