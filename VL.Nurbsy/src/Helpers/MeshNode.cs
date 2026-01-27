using System.Reactive.Disposables;
using Stride.Engine;
using Stride.Rendering;
using Stride.Rendering.ProceduralModels;
using VL.Core;
using VL.Core.Import;
using VL.Lib.Basics.Resources;
using VL.Stride;
using StrideModel = Stride.Rendering.Model;

namespace VL.Nurbsy.Helpers
{
    /// <summary>
    /// Base class for procedural mesh generation nodes.
    /// </summary>
    [ProcessNode]
    public abstract class MeshNode : IDisposable
    {
        public Mesh Output { get; protected set; }

        private readonly SerialDisposable _meshDisposable = new();

        protected void Generate()
        {
            var gameProvider = AppHost.Current.Services.GetGameProvider();
            var key = GetResourceKey();

            var provider = ResourceProvider.NewPooledSystemWide(
                key,
                _ =>
                {
                    return gameProvider.Bind(game =>
                    {
                        var model = new StrideModel();
                        var generator = Build();

                        generator.Generate(game.Services, model);

                        return ResourceProvider.Return(model.Meshes[0], DisposeMeshBuffers);
                    });
                }
            );

            var meshHandle = provider.GetHandle();
            _meshDisposable.Disposable = meshHandle;
            Output = meshHandle.Resource;
        }

        protected abstract object GetResourceKey();

        protected abstract IProceduralModel Build();

        protected static IResourceProvider<Game> GetGameProvider()
        {
            return AppHost.Current.Services.GetGameProvider();
        }

        private static void DisposeMeshBuffers(Mesh m)
        {
            if (m.Draw is null)
                return;
            m.Draw.IndexBuffer?.Buffer?.Dispose();
            foreach (var v in m.Draw.VertexBuffers)
                v.Buffer?.Dispose();
        }

        public void Dispose()
        {
            _meshDisposable.Dispose();
        }
    }
}
