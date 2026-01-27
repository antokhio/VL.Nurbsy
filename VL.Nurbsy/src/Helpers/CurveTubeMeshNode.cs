using VL.Core.Import;

namespace VL.Nurbsy.Helpers
{
    /// <summary>
    /// Base class for curve tube mesh generation.
    /// Currently: uniform tesselation.
    /// </summary>
    [ProcessNode]
    public abstract class CurveTubeMeshNode<TCurve> : MeshNode
    {
        protected const float DefaultRadius = 0.1f;
        protected const int DefaultTesselation = 32;
        protected const int DefultSegments = 12;

        protected TCurve Curve;
        protected float Radius = DefaultRadius;
        protected int Tesselation = DefaultTesselation;
        protected int Segments = DefultSegments;

        protected CurveTubeMeshNode(TCurve curve)
        {
            Curve = curve;
        }

        public virtual void SetCurve(TCurve curve)
        {
            if (!EqualityComparer<TCurve>.Default.Equals(Curve, curve))
            {
                Curve = curve;

                Generate();
            }
        }

        public void SetRadius(float radius = DefaultRadius)
        {
            if (Radius != radius)
            {
                Radius = radius;
                Generate();
            }
        }

        public void SetTessellation(int tessellation = DefaultTesselation)
        {
            if (Tesselation != tessellation)
            {
                Tesselation = tessellation;
                Generate();
            }
        }

        public void SetSegments(int segments = DefultSegments)
        {
            if (Segments != segments)
            {
                Segments = segments;
                Generate();
            }
        }
    }
}
