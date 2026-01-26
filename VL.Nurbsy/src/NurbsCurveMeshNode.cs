using System.Reactive.Disposables;
using Nurbsy;
using Nurbsy.Rendering;
using Stride.Rendering;
using VL.Core;
using VL.Core.Import;
using VL.Lib.Basics.Resources;
using VL.Nurbsy.Rendering;
using VL.Stride;
using StrideModel = Stride.Rendering.Model;

namespace VL.Nurbsy
{
    [ProcessNode]
    public abstract class NurbsCurveMeshNode<T> : IDisposable
        where T : struct
    {
        protected const float DefaultRadius = 0.1f;
        protected const int DefaultTesselation = 32;
        protected const int DeafultSegments = 12;

        private NurbsCurve<T> _curve;
        private float _radius = DefaultRadius;
        private int _tesselation = DefaultTesselation;
        private int _segments = DeafultSegments;

        public Mesh Output { get; protected set; }

        private readonly SerialDisposable _meshDisposable = new();

        protected NurbsCurveMeshNode(NurbsCurve<T> curve)
        {
            _curve = curve;
            Generate();
        }

        protected virtual void Generate()
        {
            var gameProvider = AppHost.Current.Services.GetGameProvider();

            var curve = _curve;
            var radius = _radius;
            var tessellation = _tesselation;
            var segements = _segments;

            var key = (
                gameProvider,
                typeof(BezierSurfaceModel<T>),
                curve,
                radius,
                tessellation,
                segements
            );
            var provider = ResourceProvider.NewPooledSystemWide(
                key,
                _ =>
                {
                    return gameProvider.Bind(game =>
                    {
                        var model = new StrideModel();
                        var generator = new NurbsCurveModel<T>(
                            curve,
                            radius,
                            tessellation,
                            segements
                        );

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

        public virtual void SetCurve(NurbsCurve<T> curve)
        {
            if (_curve != curve)
            {
                _curve = curve;

                Generate();
            }
        }

        public void SetRadius(float radius = DefaultRadius)
        {
            if (_radius != radius)
            {
                _radius = radius;

                Generate();
            }
        }

        public void SetTessellation(int tessellation = DefaultTesselation)
        {
            if (_tesselation != tessellation)
            {
                _tesselation = tessellation;

                Generate();
            }
        }

        public void SetSegments(int segments = DeafultSegments)
        {
            if (_segments != segments)
            {
                _segments = segments;

                Generate();
            }
        }

        public void Dispose()
        {
            _meshDisposable.Dispose();
        }
    }
}
