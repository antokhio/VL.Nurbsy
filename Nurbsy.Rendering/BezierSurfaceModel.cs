using System.Runtime.CompilerServices;
using Nurbsy;
using Nurbsy.Rendering.Helpers;
using Stride.Core.Mathematics;
using Stride.Graphics;
using Stride.Rendering.ProceduralModels;

namespace VL.Nurbsy.Rendering
{
    /// <summary>
    /// Base class to draw bezier surface via stride procedural model.
    /// </summary>
    public class BezierSurfaceModel<T> : PrimitiveProceduralModelBase
        where T : struct
    {
        public BezierSurface<T> Surface { get; set; }
        public Int2 Tesselation { get; set; }

        public BezierSurfaceModel(BezierSurface<T> surface, Int2 tesselation)
        {
            Surface = surface;
            Tesselation = tesselation;
        }

        protected override GeometricMeshData<VertexPositionNormalTexture> CreatePrimitiveMeshData()
        {
            if (typeof(T) == typeof(Vector2))
            {
                var surfaceT = Surface;
                var surface = Unsafe.As<BezierSurface<T>, BezierSurface<Vector2>>(ref surfaceT);
                return BezierSurfaceModelFactory.GenerateMeshData(surface, Tesselation);
            }

            if (typeof(T) == typeof(Vector3))
            {
                var surfaceT = Surface;
                var surface = Unsafe.As<BezierSurface<T>, BezierSurface<Vector3>>(ref surfaceT);
                return BezierSurfaceModelFactory.GenerateMeshData(surface, Tesselation);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
