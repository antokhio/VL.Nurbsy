using System.Runtime.CompilerServices;
using Nurbsy.Rendering.Helpers;
using Stride.Core.Mathematics;
using Stride.Graphics;
using Stride.Rendering.ProceduralModels;

namespace Nurbsy.Rendering
{
    public class NurbsCurveTubeModel<T> : PrimitiveProceduralModelBase
        where T : struct
    {
        public NurbsCurve<T> Curve { get; set; }
        public float Radius { get; set; }
        public int Tesslation { get; set; }
        public int Segemnts { get; set; }
        public TubeCapping Capping { get; set; }

        public NurbsCurveTubeModel(
            NurbsCurve<T> nurbsCurve,
            float radius,
            int tesselation,
            int segments
        )
        {
            Curve = nurbsCurve;
            Radius = radius;
            Tesslation = tesselation;
            Segemnts = segments;
        }

        protected override GeometricMeshData<VertexPositionNormalTexture> CreatePrimitiveMeshData()
        {
            if (typeof(T) == typeof(Vector2))
            {
                //var curveT = Curve;
                //var curve = Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref curveT);
                //return NurbsCurveModelFactory.GenerateMeshData(curve, Tesslation);
            }

            if (typeof(T) == typeof(Vector3))
            {
                var curveT = Curve;
                var curve = Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref curveT);
                return NurbsCurveTubeModelFactory.GenerateCappedTube(
                    curve,
                    Radius,
                    Tesslation,
                    Segemnts
                );
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
