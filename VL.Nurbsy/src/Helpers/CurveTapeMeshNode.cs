using System.ComponentModel;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy.Helpers
{
    [ProcessNode]
    public abstract class CurveTapeMeshNode<TCurve> : MeshNode
    {
        protected const float DefaultWidth = 0.1f;
        protected const int DefaultTesselation = 32;
        protected static readonly Vector3 DefaultDirection = Vector3.UnitZ;
        protected const string DefaultDirectionValue = "0, 0, 1";

        protected TCurve Curve;
        protected float Width = DefaultWidth;
        protected int Tesselation = DefaultTesselation;
        protected Vector3 Direction = DefaultDirection;

        protected CurveTapeMeshNode(TCurve curve)
        {
            Curve = curve;

            Generate();
        }

        public virtual void SetCurve(TCurve curve)
        {
            if (!EqualityComparer<TCurve>.Default.Equals(Curve, curve))
            {
                Curve = curve;
                Generate();
            }
        }

        public void SetWidth(float width = DefaultWidth)
        {
            if (Width != width)
            {
                Width = width;
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

        public void SetDirection([DefaultValue(DefaultDirectionValue)] Vector3 direction)
        {
            if (Direction != direction)
            {
                Direction = direction;
                Generate();
            }
        }
    }
}
