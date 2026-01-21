using Nurbsy;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// Base class for managing BezierSurface instance
    /// </summary>
    /// <typeparam name="T">Vector2, Vector3</typeparam>
    [ProcessNode(FragmentSelection = FragmentSelection.Explicit)]
    public abstract class BezierSurfaceNode<T>
        where T : struct
    {
        internal const int DefaultDegreeU = 1;
        internal const int DefaultDegreeV = 1;

        private IReadOnlyList<IReadOnlyList<ControlPoint<T>>> _controlPoints;

        private int _degreeU = DefaultDegreeU;
        private int _degreeV = DefaultDegreeV;

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

        public int DegreeU
        {
            get => _degreeU;
            protected set
            {
                _degreeU = value;
                Invalidate();
            }
        }

        public int DegreeV
        {
            get => _degreeV;
            protected set
            {
                _degreeV = value;
                Invalidate();
            }
        }

        [Fragment]
        public BezierSurface<T> Output { get; protected set; }

        public BezierSurfaceNode(BezierSurface<T> surface)
        {
            Output = surface;
        }

        public void Invalidate()
        {
            _invalidate = true;
        }

        public virtual void Build()
        {
            Output = new BezierSurface<T>(_degreeU, _degreeV, _controlPoints);
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
