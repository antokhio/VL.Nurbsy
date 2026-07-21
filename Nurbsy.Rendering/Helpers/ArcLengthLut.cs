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
    /// <see cref="NurbsSurface{T}"/>. The cache is rebuilt only when the surface's structural
    /// "shape" (degree and control point grid dimensions) changes, not when only control point
    /// positions move. This trades a small amount of arc-length accuracy under animated control
    /// points for a large reduction in per-frame CPU cost, since the expensive arc-length
    /// integration is skipped as long as the topology stays the same.
    /// </summary>
    public sealed class NurbsSurfaceArcLengthCache<T>
        where T : struct
    {
        private (int DegreeU, int DegreeV, int CountU, int CountV)? _signature;

        public int Resolution { get; set; } = ArcLengthLut.DefaultResolution;

        public ArcLengthLut ULut { get; private set; }
        public ArcLengthLut VLut { get; private set; }

        /// <summary>
        /// Returns the cached lookup tables, rebuilding them first if the surface's degree or
        /// control point grid dimensions changed since the last call.
        /// </summary>
        public void GetOrBuild(in NurbsSurface<T> surface)
        {
            int countU = surface.ControlPoints.Count;
            int countV = countU > 0 ? surface.ControlPoints[0].Count : 0;
            var signature = (surface.DegreeU, surface.DegreeV, countU, countV);

            if (ULut == null || VLut == null || _signature != signature)
            {
                Rebuild(surface);
                _signature = signature;
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
