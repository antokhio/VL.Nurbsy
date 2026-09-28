using Nurbsy;

namespace Nurbsy.Rendering.Helpers
{
    /// <summary>
    /// A precomputed lookup table mapping a normalized arc-length coordinate t in [0,1]
    /// to the real curve parameter, built by sampling <see cref="NurbsCurve{T}.GetParamAt(float)"/>
    /// at fixed intervals once and linearly interpolating between samples afterwards.
    /// This avoids repeating the (relatively expensive) Gauss-Legendre arc-length integration
    /// for every tessellation vertex.
    /// </summary>
    public sealed class ArcLengthLut
    {
        public const int DefaultResolution = 64;

        private readonly double[] _params;

        private ArcLengthLut(double[] parameters)
        {
            _params = parameters;
        }

        public static ArcLengthLut Build<T>(NurbsCurve<T> curve, int resolution)
            where T : struct
        {
            resolution = Math.Max(1, resolution);
            var parameters = new double[resolution + 1];
            for (int i = 0; i <= resolution; i++)
            {
                float t = i / (float)resolution;
                parameters[i] = curve.GetParamAt(t);
            }
            return new ArcLengthLut(parameters);
        }

        /// <summary>
        /// Resolves the real curve parameter for a normalized arc-length coordinate,
        /// linearly interpolating between the two closest precomputed samples.
        /// </summary>
        public double Sample(double t)
        {
            t = Math.Clamp(t, 0.0, 1.0);

            double scaled = t * (_params.Length - 1);
            int i0 = (int)scaled;
            if (i0 >= _params.Length - 1)
                return _params[_params.Length - 1];

            int i1 = i0 + 1;
            double frac = scaled - i0;
            return _params[i0] + (_params[i1] - _params[i0]) * frac;
        }
    }

    /// <summary>
    /// Caches the arc-length lookup tables used to correct tessellation compression for a
    /// <see cref="NurbsSurface{T}"/>. The cache is rebuilt whenever the control net, the knot
    /// vectors, the degrees or the resolution change. Change detection uses reference equality on
    /// the (immutable) control point and knot collections, so an unchanged surface reuses the
    /// tables while any edit produces exactly the same result as a freshly built cache.
    /// </summary>
    public sealed class NurbsSurfaceArcLengthCache<T>
        where T : struct
    {
        private object? _controlPoints;
        private object? _knotsU;
        private object? _knotsV;
        private (int DegreeU, int DegreeV, int Resolution)? _signature;

        public int Resolution { get; set; } = ArcLengthLut.DefaultResolution;

        public ArcLengthLut ULut { get; private set; }
        public ArcLengthLut VLut { get; private set; }

        /// <summary>
        /// Returns the cached lookup tables, rebuilding them first if the surface's control net,
        /// knots, degrees or the resolution changed since the last call.
        /// </summary>
        public void GetOrBuild(in NurbsSurface<T> surface)
        {
            var signature = (surface.DegreeU, surface.DegreeV, Resolution);

            if (ULut == null || VLut == null || _signature != signature
                || !ReferenceEquals(_controlPoints, surface.ControlPoints)
                || !ReferenceEquals(_knotsU, surface.KnotsU)
                || !ReferenceEquals(_knotsV, surface.KnotsV))
            {
                Rebuild(surface);
                _signature = signature;
                _controlPoints = surface.ControlPoints;
                _knotsU = surface.KnotsU;
                _knotsV = surface.KnotsV;
            }
        }

        private void Rebuild(in NurbsSurface<T> surface)
        {
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;

            // Reference iso-curve along U, taken at the middle of the V domain, used to map
            // normalized tu -> real u. Mirrors NurbsSurface<T>.Sample's arc-length correction.
            float vRef = (float)((knotsV[0] + knotsV[knotsV.Count - 1]) * 0.5);
            var refCurveU = surface.GetCurveV(vRef);
            ULut = ArcLengthLut.Build(refCurveU, Resolution);

            // Reference iso-curve along V, taken at the middle of the U domain, used to map
            // normalized tv -> real v. Unlike the exact per-column approach, this single reference
            // curve is reused for every column, which is the approximation that makes caching possible.
            float uRef = (float)((knotsU[0] + knotsU[knotsU.Count - 1]) * 0.5);
            var refCurveV = surface.GetCurveU(uRef);
            VLut = ArcLengthLut.Build(refCurveV, Resolution);
        }
    }
}
