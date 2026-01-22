using Nurbsy;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// Base class for managing <see cref="BezierCurve{T}"/> instance
    /// </summary>
    /// <typeparam name="T">Vector2, Vector3</typeparam>
    [ProcessNode(FragmentSelection = FragmentSelection.Explicit)]
    public abstract class BezierCurveNode<T>
        where T : struct
    {
        private IReadOnlyList<ControlPoint<T>> _controlPoints;

        private int? _degree;

        private bool _invalidate = false;

        public IReadOnlyList<ControlPoint<T>> ControlPoints
        {
            get => _controlPoints;
            protected set
            {
                _controlPoints = value;
                Invalidate();
            }
        }

        public int? Degree
        {
            get => _degree ?? Math.Max(1, _controlPoints.Count - 1);
            protected set
            {
                _degree = value;
                Invalidate();
            }
        }

        [Fragment]
        public BezierCurve<T> Output { get; protected set; }

        protected BezierCurveNode(BezierCurve<T> curve)
        {
            // We populate only control poins here since
            // they miss default value in abstract class
            _controlPoints = curve.ControlPoints;

            Output = curve;
        }

        public void Invalidate()
        {
            _invalidate = true;
        }

        public virtual void Build()
        {
            // Resolve effective degrees
            var degree = _degree ?? Math.Max(1, _controlPoints.Count - 1);

            Output = new BezierCurve<T>(degree, _controlPoints);
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
