using System.Runtime.CompilerServices;
using Nurbsy.Helpers;
using Stride.Core.Mathematics;

namespace Nurbsy
{
    public record struct NurbsSurface<T>
        where T : struct
    {
        public int DegreeU { get; }
        public int DegreeV { get; }

        public IReadOnlyList<IReadOnlyList<ControlPoint<T>>> ControlPoints { get; }
        public IReadOnlyList<double> KnotsU { get; }
        public IReadOnlyList<double> KnotsV { get; }

        public NurbsSurface(
            int degreeU,
            int degreeV,
            IReadOnlyList<IReadOnlyList<ControlPoint<T>>> controlPoints,
            IReadOnlyList<double> knotsU,
            IReadOnlyList<double> knotsV
        )
        {
            DegreeU = degreeU;
            DegreeV = degreeV;
            ControlPoints = controlPoints;
            KnotsU = knotsU;
            KnotsV = knotsV;

            Check();
        }

        public void Check()
        {
            Validate.Argument(DegreeU > 0, nameof(DegreeU), "Degree must be greater than zero.");
            Validate.Argument(DegreeV > 0, nameof(DegreeV), "Degree must be greater than zero.");

            Validate.Argument(
                KnotsU.Count > 0,
                nameof(KnotsU),
                "KnotVector size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(KnotsU),
                nameof(KnotsU),
                "KnotVector must be a non-decreasing sequence of real numbers."
            );

            Validate.Argument(
                KnotsV.Count > 0,
                nameof(KnotsV),
                "KnotVector size must be greater than zero."
            );
            Validate.Argument(
                Validate.IsValidKnots(KnotsV),
                nameof(KnotsV),
                "KnotVector must be a non-decreasing sequence of real numbers."
            );

            Validate.Argument(
                ControlPoints.Count > 0,
                nameof(ControlPoints),
                "ControlPoints must contain one point at least."
            );
            Validate.Argument(
                Validate.IsValidNURBS(DegreeU, ControlPoints.Count, KnotsU.Count),
                nameof(ControlPoints),
                "Arguments must be fit: m = n + p + 1"
            );

            if (ControlPoints.Count > 0)
            {
                Validate.Argument(
                    Validate.IsValidNURBS(DegreeV, ControlPoints[0].Count, KnotsV.Count),
                    nameof(ControlPoints),
                    "Arguments must be fit: m = n + p + 1"
                );
            }
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page103
        /// Algorithm A3.5
        /// Compute surface point.
        /// </summary>
        /// <param name="uv">The uv parameter.</param>
        /// <returns>Point on surface</returns>
        public T GetPointOnSurface(Vector2 uv)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.GetPointOnSurface(in surface, uv);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.GetPointOnSurface(in surface, uv);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page111
        /// Algorithm A3.6
        /// Compute surface derivatives.
        /// </summary>
        /// <param name="derivative">The maximum derivative order to compute.</param>
        /// <param name="uv">The uv parameter.</param>
        /// <returns>2D array where result[k][l] is the mixed partial derivative ∂^(k+l)S/∂u^k∂v^l.</returns>
        public T[][] GetDerivatives(int derivative, Vector2 uv)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.ComputeDerivatives(in surface, derivative, uv);
                return Unsafe.As<Vector2[][], T[][]>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.ComputeDerivatives(in surface, derivative, uv);
                return Unsafe.As<Vector3[][], T[][]>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Optimized computation for first-order surface derivatives.
        /// </summary>
        /// <param name="uv">The uv parameter.</param>
        /// <returns>2x2 array where result[k][l] is the mixed partial derivative ∂^(k+l)S/∂u^k∂v^l for k,l ∈ {0,1}.</returns>
        public T[][] GetFirstOrderDerivatives(Vector2 uv)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.ComputeFirstOrderDerivatives(in surface, uv);
                return Unsafe.As<Vector2[][], T[][]>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.ComputeFirstOrderDerivatives(in surface, uv);
                return Unsafe.As<Vector3[][], T[][]>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }
    }
}
