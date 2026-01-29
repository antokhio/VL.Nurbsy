using Nurbsy;
using Nurbsy.Algorithm;
using Nurbsy.Helpers;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy
{
    [ProcessNode(Name = "NurbsSurface (3d Ruled)")]
    public class NurbsSurfaceNodeRuledSurface : NurbsSurfaceNode<Vector3>
    {
        protected static readonly NurbsCurve<Vector3> DefaultCurve0 = new(
            1,
            [new Vector3(-0.5f, 0.5f, 0f), new Vector3(0.5f, 0.5f, 0f)],
            KnotsUtils.GenerateClampedKnots(1, 2)
        );
        protected static readonly NurbsCurve<Vector3> DefaultCurve1 = new(
            1,
            [new Vector3(-0.5f, -0.5f, 0f), new Vector3(0.5f, -0.5f, 0f)],
            KnotsUtils.GenerateClampedKnots(1, 2)
        );

        private NurbsCurve<Vector3> _curve0;
        private NurbsCurve<Vector3> _curve1;

        public NurbsSurfaceNodeRuledSurface()
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
