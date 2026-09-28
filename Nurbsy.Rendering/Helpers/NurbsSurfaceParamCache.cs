namespace Nurbsy.Rendering.Helpers
{
    /// <summary>
    /// Caches the tessellation parameters produced by <see cref="TextureParametrization"/> for a
    /// <see cref="NurbsSurface{T}"/>. They only depend on degrees, knot vectors, control grid
    /// dimensions and tessellation, so animating control point positions reuses the cached
    /// values. Knot change detection uses reference equality on the immutable knot collections.
    /// </summary>
    public sealed class NurbsSurfaceParamCache<T>
        where T : struct
    {
        private object? _knotsU;
        private object? _knotsV;
        private (int DegreeU, int DegreeV, int CountU, int CountV, int SamplesU, int SamplesV)? _signature;

        public float[] UParams { get; private set; } = [];
        public float[] VParams { get; private set; } = [];

        /// <summary>
        /// Ensures <see cref="UParams"/> and <see cref="VParams"/> match the given surface and
        /// sample counts, rebuilding them only when the surface topology changed.
        /// </summary>
        public void GetOrBuild(in NurbsSurface<T> surface, int samplesU, int samplesV)
        {
            var controlPoints = surface.ControlPoints;
            var signature = (
                surface.DegreeU,
                surface.DegreeV,
                controlPoints.Count,
                controlPoints[0].Count,
                samplesU,
                samplesV
            );

            if (_signature != signature
                || !ReferenceEquals(_knotsU, surface.KnotsU)
                || !ReferenceEquals(_knotsV, surface.KnotsV))
            {
                UParams = TextureParametrization.Resolve(
                    surface.DegreeU,
                    surface.KnotsU,
                    controlPoints.Count,
                    samplesU
                );
                VParams = TextureParametrization.Resolve(
                    surface.DegreeV,
                    surface.KnotsV,
                    controlPoints[0].Count,
                    samplesV
                );

                _signature = signature;
                _knotsU = surface.KnotsU;
                _knotsV = surface.KnotsV;
            }
        }
    }
}