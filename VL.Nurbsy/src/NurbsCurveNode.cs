using Nurbsy;
using Nurbsy.Algorithm;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// Base class for managing NurbsCurve instance
    /// </summary>
    [ProcessNode(FragmentSelection = FragmentSelection.Explicit)]
    public abstract class NurbsCurveNode<T>
        where T : struct
    {
        internal const int DefaultDegree = 1;

        private IReadOnlyList<ControlPoint<T>> _controlPoints;
        private int _degree = DefaultDegree;
        private IReadOnlyList<double> _knots;

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

        public int Degree
        {
            get => _degree;
            protected set
            {
                _degree = value;
                Invalidate();
            }
        }

        public IReadOnlyList<double> Knots
        {
            get => _knots;
            protected set
            {
                _knots = value;
                Invalidate();
            }
        }

        [Fragment]
        public NurbsCurve<T> Output { get; protected set; }

        [Fragment]
        protected NurbsCurveNode(NurbsCurve<T> curve)
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
            var knots = _knots ?? KnotsUtils.GenerateClampedKnots(_degree, _controlPoints.Count);

            Output = new NurbsCurve<T>(_degree, _controlPoints, knots);
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
