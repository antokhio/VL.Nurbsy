using Nurbsy;
using Nurbsy.Algorithm;
using VL.Core.Import;

namespace VL.Nurbsy
{
    /// <summary>
    /// Base class for managing <see cref="NurbsSurface{T}"/> instance
    /// </summary>
    /// <typeparam name="T">Vector2, Vector3</typeparam>
    [ProcessNode(FragmentSelection = FragmentSelection.Explicit)]
    public abstract class NurbsSurfaceNode<T>
        where T : struct
    {
        protected const int DefaultDegreeU = 1;
        protected const int DefaultDegreeV = 1;

        private IReadOnlyList<IReadOnlyList<ControlPoint<T>>> _controlPoints;

        private int _degreeU = DefaultDegreeU;
        private int _degreeV = DefaultDegreeV;

        private IReadOnlyList<double> _knotsU;
        private IReadOnlyList<double> _knotsV;

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

        public IReadOnlyList<double> KnotsU
        {
            get => _knotsU;
            protected set
            {
                _knotsU = value;
                Invalidate();
            }
        }

        public IReadOnlyList<double> KnotsV
        {
            get => _knotsV;
            protected set
            {
                _knotsV = value;
                Invalidate();
            }
        }

        [Fragment]
        public NurbsSurface<T> Output { get; protected set; }

        protected NurbsSurfaceNode(NurbsSurface<T> surface)
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
            var knotsU = _knotsU ?? KnotsUtils.GenerateClampedKnots(_degreeU, _controlPoints.Count);
            var knotsV =
                _knotsV ?? KnotsUtils.GenerateClampedKnots(_degreeV, _controlPoints[0].Count);

            Output = new NurbsSurface<T>(_degreeU, _degreeV, _controlPoints, knotsU, knotsV);
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
