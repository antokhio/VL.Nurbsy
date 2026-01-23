using System.ComponentModel;
using System.Reactive.Disposables;
using Nurbsy;
using Stride.Core.Mathematics;
using Stride.Rendering;
using VL.Core;
using VL.Core.Import;
using VL.Lib.Basics.Resources;
using VL.Nurbsy.Stride;
using VL.Stride;
using StrideModel = Stride.Rendering.Model;

namespace VL.Nurbsy
{
    [ProcessNode]
    public abstract class BezierSurfaceMeshNode<T> : IDisposable
        where T : struct
    {
        protected static readonly Int2 DefaultTesselation = new Int2(32, 32);
        protected const string DefaultTesselationValue = "32, 32";

        private BezierSurface<T> _surface;
        private Int2 _tesselation = DefaultTesselation;

        public Mesh Output { get; protected set; }

        private readonly SerialDisposable _meshDisposable = new();

        protected BezierSurfaceMeshNode(BezierSurface<T> surface)
        {
            _surface = surface;

            Generate();
        }

        protected virtual void Generate()
        {
            var gameProvider = AppHost.Current.Services.GetGameProvider();

            var surface = _surface;
            var tessellation = _tesselation;

            var key = (gameProvider, typeof(BezierSurfaceModel<T>), surface, tessellation);
            var provider = ResourceProvider.NewPooledSystemWide(
                key,
                _ =>
                {
                    return gameProvider.Bind(game =>
                    {
                        var model = new StrideModel();
                        var generator = new BezierSurfaceModel<T>(surface, tessellation);

                        generator.Generate(game.Services, model);

                        return ResourceProvider.Return(
                            model.Meshes[0],
                            m =>
                            {
                                // TODO: remove
                                // should be just
                                // return ResourceProvider.Return(model.Meshes[0], m => m.ReleaseGraphicsResources());

                                if (m.Draw is null)
                                    return;

                                m.Draw.IndexBuffer?.Buffer?.Dispose();
                                foreach (var v in m.Draw.VertexBuffers)
                                    v.Buffer?.Dispose();
                            }
                        );
                    });
                }
            );

            var meshHandle = provider.GetHandle();
            _meshDisposable.Disposable = meshHandle;

            Output = meshHandle.Resource;
        }

        public virtual void SetSurface(BezierSurface<T> surface)
        {
            if (_surface != surface)
            {
                _surface = surface;

                Generate();
            }
        }

        public void SetTessellation([DefaultValue(DefaultTesselationValue)] Int2 tessellation)
        {
            if (_tesselation != tessellation)
            {
                _tesselation = tessellation;

                Generate();
            }
        }

        public void Dispose()
        {
            _meshDisposable.Dispose();
        }
    }
}
