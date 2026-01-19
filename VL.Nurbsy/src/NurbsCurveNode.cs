using Nurbsy;
using Stride.Core.Mathematics;
using VL.Core.Import;
using VL.Model;

namespace VL.Nurbsy
{
    [ProcessNode]
    public abstract class NurbsCurveNode<T>
        where T : struct
    {
        internal const int DefaultDegree = 2;

        private NurbsCurve<T> _output;
        protected bool _invalidated = false;

        protected int _degree = DefaultDegree;
        protected IReadOnlyList<T> _controlPoints;
        protected IReadOnlyList<float> _knots;

        public NurbsCurve<T> Output => _output;

        protected NurbsCurveNode(
            int defaultDegree,
            IReadOnlyList<T> defaultControlPoints,
            IReadOnlyList<float> defaultKnots = null
        )
        {
            _degree = defaultDegree;
            _controlPoints = defaultControlPoints;
            _knots = defaultKnots;

            // Initialize immediately
            Rebuild();
        }

        public void SetDegree(int degree = DefaultDegree)
        {
            if (_degree != degree)
            {
                _degree = degree;
                _invalidated = true;
            }
        }

        public void SetControlPoints(IReadOnlyList<T> controlPoints)
        {
            if (_controlPoints != controlPoints)
            {
                _controlPoints = controlPoints;
                _invalidated = true;
            }
        }

        public void SetKnots([Pin(Visibility = PinVisibility.Optional)] IReadOnlyList<float> knots)
        {
            if (_knots != knots)
            {
                _knots = knots;
                _invalidated = true;
            }
        }

        public void Update()
        {
            if (_invalidated)
            {
                Rebuild();
                _invalidated = false;
            }
        }

        private void Rebuild()
        {
            try
            {
                if (_knots != null)
                {
                    var knotsDouble = _knots.Select(k => (double)k).ToArray();
                    _output = new NurbsCurve<T>(_degree, _controlPoints, knotsDouble);
                }
                else
                {
                    var controlPoints = _controlPoints
                        .Select(p => new ControlPoint<T>(p))
                        .ToArray();
                    _output = new NurbsCurve<T>(_degree, controlPoints);
                }
            }
            catch
            {
                // Fallback or handle error (e.g. invalid arguments)
                _output = default;
                throw;
            }
        }
    }

    [ProcessNode(Name = "NurbsCurve (2D)")]
    public class NurbsCurve2DNode : NurbsCurveNode<Vector2>
    {
        static readonly IReadOnlyList<Vector2> DefaultControlPoints =
        [
            new(-0.5f, 0f),
            new(0f, 0f),
            new(0.5f, 0f),
        ];

        public NurbsCurve2DNode()
            : base(DefaultDegree, DefaultControlPoints) { }
    }

    [ProcessNode(Name = "NurbsCurve (3D)")]
    public class NurbsCurve3DNode : NurbsCurveNode<Vector3>
    {
        static readonly IReadOnlyList<Vector3> DefaultControlPoints =
        [
            new(-0.5f, 0f, 0f),
            new(0f, 0f, 0f),
            new(0.5f, 0f, 0f),
        ];

        public NurbsCurve3DNode()
            : base(DefaultDegree, DefaultControlPoints) { }
    }
}
