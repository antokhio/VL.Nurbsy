using Nurbsy;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// Base class for managing <see cref="BezierSurface{T}"/> instance
    /// </summary>
    /// <typeparam name="T">Vector2, Vector3</typeparam>
    [ProcessNode]
    public abstract class BezierSurfaceNode<T>
        where T : struct
    {
        private IReadOnlyList<IReadOnlyList<ControlPoint<T>>> _controlPoints;

        // We calculate degree from control points if not set
        private int? _degreeU;
        private int? _degreeV;

        private bool _invalidate = false;

        public IReadOnlyList<IReadOnlyList<ControlPoint<T>>> ControlPoints
        {
            get => _controlPoints;
            protected set
            {
                _controlPoints = value;
                Invalidate();
            }
        }

        public int? DegreeU
        {
            get => _degreeU ?? Math.Max(1, _controlPoints.Count - 1);
            protected set
            {
                _degreeU = value;
                Invalidate();
            }
        }
        public int? DegreeV
        {
            get => _degreeV ?? Math.Max(1, _controlPoints[0].Count - 1);
            protected set
            {
                _degreeV = value;
                Invalidate();
            }
        }

        [Fragment]
        public BezierSurface<T> Output { get; protected set; }

        protected BezierSurfaceNode(BezierSurface<T> surface)
        {
            // We populate only control poins here since
            // they miss default value in abstract class
            _controlPoints = surface.ControlPoints;

            Output = surface;
        }

        public void Invalidate()
        {
            _invalidate = true;
        }

        public virtual void Build()
        {
            // Resolve effective degrees
            var degreeU = _degreeU ?? Math.Max(1, _controlPoints.Count - 1);
            var degreeV = _degreeV ?? Math.Max(1, _controlPoints[0].Count - 1);

            Output = new BezierSurface<T>(degreeU, degreeV, _controlPoints);
        }

        [Fragment(Order = int.MaxValue)]
        public void Update()
        {
            if (_invalidate)
            {
                Build();
                _invalidate = false;
            }
        }
    }
}
