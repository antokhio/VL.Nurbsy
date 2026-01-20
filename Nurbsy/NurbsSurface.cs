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

        /// <summary>
        /// The NURBS Book 2nd Edition Page114
        /// Algorithm A3.7
        /// Compute control points of derivative surfaces.
        /// </summary>
        /// <param name="derivative">The maximum derivative order.</param>
        /// <param name="minSpanIndexU">Minimum span index in U direction.</param>
        /// <param name="maxSpanIndexU">Maximum span index in U direction.</param>
        /// <param name="minSpanIndexV">Minimum span index in V direction.</param>
        /// <param name="maxSpanIndexV">Maximum span index in V direction.</param>
        /// <returns>4D array PKL[k][l][i][j] of derivative control points.</returns>
        public T[][][][] GetControlPointsOfDerivatives(
            int derivative,
            int minSpanIndexU,
            int maxSpanIndexU,
            int minSpanIndexV,
            int maxSpanIndexV
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.ComputeControlPointsOfDerivatives(
                    in surface,
                    derivative,
                    minSpanIndexU,
                    maxSpanIndexU,
                    minSpanIndexV,
                    maxSpanIndexV
                );
                return Unsafe.As<Vector2[][][][], T[][][][]>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.ComputeControlPointsOfDerivatives(
                    in surface,
                    derivative,
                    minSpanIndexU,
                    maxSpanIndexU,
                    minSpanIndexV,
                    maxSpanIndexV
                );
                return Unsafe.As<Vector3[][][][], T[][][][]>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page115
        /// Algorithm A3.8
        /// Compute surface derivatives using all basis functions.
        /// </summary>
        /// <param name="derivative">The maximum derivative order to compute.</param>
        /// <param name="uv">The uv parameter.</param>
        /// <returns>2D array where result[k][l] is the mixed partial derivative ∂^(k+l)S/∂u^k∂v^l.</returns>
        public T[][] GetDerivativesByAllBasisFunctions(int derivative, Vector2 uv)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.ComputeDerivativesByAllBasisFunctions(
                    in surface,
                    derivative,
                    uv
                );
                return Unsafe.As<Vector2[][], T[][]>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.ComputeDerivativesByAllBasisFunctions(
                    in surface,
                    derivative,
                    uv
                );
                return Unsafe.As<Vector3[][], T[][]>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 137, Algorithm A4.4.
        /// Compute rational surface derivatives S(u,v).
        /// Properly handles the quotient rule for NURBS (rational) surfaces.
        /// </summary>
        /// <param name="derivative">The maximum derivative order to compute.</param>
        /// <param name="uv">The uv parameter.</param>
        /// <returns>2D array where result[k][l] is the mixed partial derivative ∂^(k+l)S/∂u^k∂v^l.</returns>
        public T[][] GetRationalDerivatives(int derivative, Vector2 uv)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.ComputeRationalSurfaceDerivatives(
                    in surface,
                    derivative,
                    uv
                );
                return Unsafe.As<Vector2[][], T[][]>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.ComputeRationalSurfaceDerivatives(
                    in surface,
                    derivative,
                    uv
                );
                return Unsafe.As<Vector3[][], T[][]>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Optimized computation of rational surface first-order derivatives.
        /// Returns S(u,v), Su, and Sv using the quotient rule.
        /// </summary>
        /// <param name="uv">The uv parameter.</param>
        /// <param name="S">Output: Surface point S(u,v).</param>
        /// <param name="Su">Output: First derivative in U direction ∂S/∂u.</param>
        /// <param name="Sv">Output: First derivative in V direction ∂S/∂v.</param>
        public void GetRationalFirstOrderDerivatives(Vector2 uv, out T S, out T Su, out T Sv)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                NurbsSurface2DHelper.ComputeRationalSurfaceFirstOrderDerivatives(
                    in surface,
                    uv,
                    out var s2,
                    out var su2,
                    out var sv2
                );
                S = Unsafe.As<Vector2, T>(ref s2);
                Su = Unsafe.As<Vector2, T>(ref su2);
                Sv = Unsafe.As<Vector2, T>(ref sv2);
                return;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                NurbsSurface3DHelper.ComputeRationalSurfaceFirstOrderDerivatives(
                    in surface,
                    uv,
                    out var s3,
                    out var su3,
                    out var sv3
                );
                S = Unsafe.As<Vector3, T>(ref s3);
                Su = Unsafe.As<Vector3, T>(ref su3);
                Sv = Unsafe.As<Vector3, T>(ref sv3);
                return;
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Compute surface normal at the given UV parameter.
        /// For 2D surfaces, returns zero vector.
        /// </summary>
        /// <param name="uv">The UV parameter.</param>
        /// <returns>The normalized surface normal vector.</returns>
        public T GetNormal(Vector2 uv)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.ComputeNormal(in surface, uv);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.ComputeNormal(in surface, uv);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Compute surface curvature at the given UV parameter.
        /// For 2D surfaces, returns 0.
        /// </summary>
        /// <param name="curvature">The type of curvature to compute.</param>
        /// <param name="uv">The UV parameter.</param>
        /// <returns>The curvature value.</returns>
        public double GetCurvature(SurfaceCurvature curvature, Vector2 uv)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                return NurbsSurface2DHelper.ComputeCurvature(in surface, curvature, uv);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                return NurbsSurface3DHelper.ComputeCurvature(in surface, curvature, uv);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Swap the U and V directions of the surface.
        /// Creates a new surface where U becomes V and V becomes U.
        /// </summary>
        /// <returns>A new surface with swapped U and V directions.</returns>
        public NurbsSurface<T> SwapDim()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.SwapDim(in surface);
                return Unsafe.As<NurbsSurface<Vector2>, NurbsSurface<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.SwapDim(in surface);
                return Unsafe.As<NurbsSurface<Vector3>, NurbsSurface<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Reverse the surface in the specified direction.
        /// </summary>
        /// <param name="direction">The direction to reverse (U, V, or All).</param>
        /// <returns>A new reversed surface.</returns>
        public NurbsSurface<T> Reverse(SurfaceDirection direction = SurfaceDirection.All)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.Reverse(in surface, direction);
                return Unsafe.As<NurbsSurface<Vector2>, NurbsSurface<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.Reverse(in surface, direction);
                return Unsafe.As<NurbsSurface<Vector3>, NurbsSurface<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 137, Algorithm A5.3.
        /// Insert a knot into the surface along U or V direction.
        /// </summary>
        /// <param name="insertKnot">The knot value to insert.</param>
        /// <param name="times">Number of times to insert the knot.</param>
        /// <param name="direction">The direction to insert (UDirection or VDirection).</param>
        /// <param name="result">The resulting surface with inserted knot.</param>
        /// <returns>The number of times the knot was actually inserted.</returns>
        public int InsertKnot(
            double insertKnot,
            int times,
            SurfaceDirection direction,
            out NurbsSurface<T> result
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                int inserted = NurbsSurface2DHelper.InsertKnot(
                    in surface,
                    insertKnot,
                    times,
                    direction,
                    out var result2
                );
                result = Unsafe.As<NurbsSurface<Vector2>, NurbsSurface<T>>(ref result2);
                return inserted;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                int inserted = NurbsSurface3DHelper.InsertKnot(
                    in surface,
                    insertKnot,
                    times,
                    direction,
                    out var result3
                );
                result = Unsafe.As<NurbsSurface<Vector3>, NurbsSurface<T>>(ref result3);
                return inserted;
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Refine the knot vector by inserting multiple knots.
        /// </summary>
        /// <param name="insertKnotElements">The knots to insert.</param>
        /// <param name="direction">The direction to refine (UDirection or VDirection).</param>
        /// <returns>A new surface with refined knot vector.</returns>
        public NurbsSurface<T> RefineKnotVector(
            IReadOnlyList<double> insertKnotElements,
            SurfaceDirection direction
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.RefineKnotVector(
                    in surface,
                    insertKnotElements,
                    direction
                );
                return Unsafe.As<NurbsSurface<Vector2>, NurbsSurface<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.RefineKnotVector(
                    in surface,
                    insertKnotElements,
                    direction
                );
                return Unsafe.As<NurbsSurface<Vector3>, NurbsSurface<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }
    }
}
