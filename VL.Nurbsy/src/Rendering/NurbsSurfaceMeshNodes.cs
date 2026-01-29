using Nurbsy;
using Nurbsy.Rendering;
using Stride.Core.Mathematics;
using Stride.Rendering.ProceduralModels;
using VL.Core.Import;
using VL.Nurbsy.Helpers;

namespace VL.Nurbsy.Rendering
{
    /// <summary>
    /// 2D Nurbs surface mesh generator.
    /// Excludes normals computation.
    /// </summary>
    [ProcessNode(Name = "NurbsSurface (2d Mesh)")]
    public class NurbsSurfaceMeshNode2D : SurfaceMeshNode<NurbsSurface<Vector2>>
    {
        public NurbsSurfaceMeshNode2D()
            : base(NurbsSurfaceNode2D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsSurfaceModel<Vector2>(Surface, Tesselation);
        }

        protected override object GetResourceKey()
        {
            return (GetGameProvider(), typeof(NurbsSurfaceModel<Vector2>), Surface, Tesselation);
        }
    }

    /// <summary>
    /// 3D Bezier surface mesh generator.
    /// Includes normals computation.
    /// </summary>
    [ProcessNode(Name = "NurbsSurface (3d Mesh)")]
    public class NurbsSurfaceMeshNode3D : SurfaceMeshNode<NurbsSurface<Vector3>>
    {
        public NurbsSurfaceMeshNode3D()
            : base(NurbsSurfaceNode3D.CreateDefault()) { }

        protected override IProceduralModel Build()
        {
            return new NurbsSurfaceModel<Vector3>(Surface, Tesselation);
        }

        protected override object GetResourceKey()
        {
            return (GetGameProvider(), typeof(NurbsSurfaceModel<Vector3>), Surface, Tesselation);
        }
    }
}
