using System.Runtime.CompilerServices;
using Nurbsy.Rendering.Helpers;
using Stride.Core.Mathematics;
using Stride.Graphics;
using Stride.Rendering.ProceduralModels;

namespace Nurbsy.Rendering
{
    /// <summary>
    /// Base class for NURBS surface model generation.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class NurbsSurfaceModel<T> : PrimitiveProceduralModelBase
        where T : struct
    {
        public NurbsSurface<T> Surface { get; set; }

        public Int2 Tesselation { get; set; }

        public NurbsSurfaceModel(NurbsSurface<T> surface, Int2 tesselation)
        {
            Surface = surface;
            Tesselation = tesselation;
        }

        protected override GeometricMeshData<VertexPositionNormalTexture> CreatePrimitiveMeshData()
        {
            if (typeof(T) == typeof(Vector2))
            {
                var surfaceT = Surface;
                var surface = Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref surfaceT);
                return NurbsSurfaceModelFactory.GenerateMeshData(surface, Tesselation);
            }

            if (typeof(T) == typeof(Vector3))
            {
                var surfaceT = Surface;
                var surface = Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref surfaceT);
                return NurbsSurfaceModelFactory.GenerateMeshData(surface, Tesselation);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
