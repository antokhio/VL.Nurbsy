using Nurbsy;
using Nurbsy.Algorithm;
using Nurbsy.Helpers;
using Stride.Core.Mathematics;
using VL.Core;
using VL.Core.Import;
using VL.Lib.Collections;

namespace VL.Nurbsy
{
    [ProcessNode(Name = "NurbsSurface (3d Loft)")]
    public class NurbsSurfaceNodeLoft : NurbsSurfaceNode<Vector3>
    {
        protected static readonly IReadOnlyList<NurbsCurve<Vector3>> DefaultSections =
        [
            new(
                1,
                [new Vector3(-0.5f, 0.5f, 0f), new Vector3(-0.5f, 0.5f, 0f)],
                KnotsUtils.GenerateClampedKnots(1, 2)
            ),
            new(
                1,
                [new Vector3(-0.5f, -0.5f, 0f), new Vector3(-0.5f, -0.5f, 0f)],
                KnotsUtils.GenerateClampedKnots(1, 2)
            ),
        ];

        private IReadOnlyList<NurbsCurve<Vector3>> _sections;
        private Optional<int> _degreeV;
        private Optional<IReadOnlyList<double>> customTrajectoryParameters;

        public NurbsSurfaceNodeLoft()
            : base(NurbsSurfaceHelper.CreateLoftSurface(DefaultSections)) { }

        public virtual void SetSections(Spread<NurbsCurve<Vector3>> sections)
        {
            if (_sections != sections)
            {
                _sections = sections;
                Invalidate();
            }
        }

        public void SetDegreeV(Optional<int> degreeV)
        {
            if (_degreeV != degreeV)
            {
                _degreeV = degreeV;
                Invalidate();
            }
        }

        public void SetCustomTrajectoryParameters(Optional<IReadOnlyList<double>> parameters)
        {
            if (customTrajectoryParameters != parameters)
            {
                customTrajectoryParameters = parameters;
                Invalidate();
            }
        }

        protected override void Build()
        {
            Output = NurbsSurfaceHelper.CreateLoftSurface(
                _sections,
                _degreeV.HasValue ? _degreeV.Value : null,
                customTrajectoryParameters.HasValue ? customTrajectoryParameters.Value : null
            );
        }
    }
}
