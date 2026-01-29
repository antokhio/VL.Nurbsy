using System.Runtime.CompilerServices;
using Nurbsy.Rendering.Helpers;
using Stride.Core.Mathematics;
using Stride.Graphics;
using Stride.Rendering.ProceduralModels;

namespace Nurbsy.Rendering
{
    public class NurbsCurveTapeModel<T> : PrimitiveProceduralModelBase
        where T : struct
    {
        public NurbsCurve<T> Curve { get; set; }
        public float Width { get; set; }
        public int Tesselation { get; set; }
        public Vector3 Direction { get; set; } = Vector3.UnitZ;

        public NurbsCurveTapeModel(
            NurbsCurve<T> nurbsCurve,
            float width,
            int tesselation,
            Vector3? direction
        )
        {
            Curve = nurbsCurve;
            Width = width;
            Tesselation = tesselation;
            Direction = direction ?? Direction;
        }

        protected override GeometricMeshData<VertexPositionNormalTexture> CreatePrimitiveMeshData()
        {
            if (typeof(T) == typeof(Vector2))
            {
                var curveT = Curve;
                var curve = Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref curveT);
                return NurbsCurveTapeModelFactory.GenerateTape(curve, Width, Tesselation);
            }
            if (typeof(T) == typeof(Vector3))
            {
                var curveT = Curve;
                var curve = Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref curveT);
                return NurbsCurveTapeModelFactory.GenerateTape(
                    curve,
                    Width,
                    Tesselation,
                    Direction
                );
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
