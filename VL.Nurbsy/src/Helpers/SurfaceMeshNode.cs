using System.ComponentModel;
using Stride.Core.Mathematics;
using VL.Core.Import;

namespace VL.Nurbsy.Helpers
{
    /// <summary>
    /// Base class for surface generation nodes.
    /// </summary>
    [ProcessNode]
    public abstract class SurfaceMeshNode<TSurface> : MeshNode
    {
        protected static readonly Int2 DefaultTesselation = new Int2(32, 32);
        protected const string DefaultTesselationValue = "32, 32";

        protected TSurface Surface;
        protected Int2 Tesselation = DefaultTesselation;

        protected SurfaceMeshNode(TSurface surface)
        {
            Surface = surface;
        }

        public virtual void SetSurface(TSurface surface)
        {
            if (!EqualityComparer<TSurface>.Default.Equals(Surface, surface))
            {
                Surface = surface;

                Generate();
            }
        }

        public void SetTessellation([DefaultValue(DefaultTesselationValue)] Int2 tessellation)
        {
            if (Tesselation != tessellation)
            {
                Tesselation = tessellation;

                Generate();
            }
        }
    }
}
