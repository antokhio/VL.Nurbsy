/*
 * Ported by:
 * 2025 - Anton Kalabuhov (antokhio)
 *
 * Original Author:
 * 2023 - Yuqing Liang (BIMCoder Liang)
 *
 * Based on LNLib: https://github.com/BIMCoderLiang/LNLib
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using System.Collections.Immutable;
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

        public NurbsSurface(
            int degreeU,
            int degreeV,
            IReadOnlyList<IReadOnlyList<T>> controlPoints,
            IReadOnlyList<double> knotsU,
            IReadOnlyList<double> knotsV
        )
        {
            DegreeU = degreeU;
            DegreeV = degreeV;
            ControlPoints = controlPoints
                .Select(row =>
                    (IReadOnlyList<ControlPoint<T>>)
                        row.Select(cp => new ControlPoint<T>(cp)).ToImmutableArray()
                )
                .ToImmutableArray();
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
                nameof(KnotsU),
                "Arguments must be fit: m = n + p + 1"
            );

            if (ControlPoints.Count > 0)
            {
                Validate.Argument(
                    Validate.IsValidNURBS(DegreeV, ControlPoints[0].Count, KnotsV.Count),
                    nameof(KnotsV),
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
        /// The NURBS Book 2nd Edition Page137
        /// Algorithm A5.3
        /// Surface knot insertion along U or V direction.
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
        /// The NURBS Book 2nd Edition Page167
        /// Algorithm A5.5
        /// Refine surface knot vector.
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

        /// <summary>
        /// The NURBS Book 2nd Edition Page177
        /// Algorithm A5.7
        /// Decompose surface into Bezier patches.
        /// </summary>
        public IReadOnlyList<BezierSurface<T>> DecomposeToBeziers()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.DecomposeToBeziers(in surface);
                return Unsafe.As<
                    IReadOnlyList<BezierSurface<Vector2>>,
                    IReadOnlyList<BezierSurface<T>>
                >(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.DecomposeToBeziers(in surface);
                return Unsafe.As<
                    IReadOnlyList<BezierSurface<Vector3>>,
                    IReadOnlyList<BezierSurface<T>>
                >(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page186
        /// Surface knot removal.
        /// </summary>
        /// <param name="removeKnot">The knot value to remove.</param>
        /// <param name="times">Number of times to remove the knot.</param>
        /// <param name="direction">The direction to remove from (UDirection or VDirection).</param>
        /// <param name="result">The resulting surface with removed knot.</param>
        /// <returns>True if the knot was successfully removed.</returns>
        public bool TryRemoveKnot(
            double removeKnot,
            int times,
            SurfaceDirection direction,
            out NurbsSurface<T> result
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                bool success = NurbsSurface2DHelper.RemoveKnot(
                    in surface,
                    removeKnot,
                    times,
                    direction,
                    out var result2
                );
                result = Unsafe.As<NurbsSurface<Vector2>, NurbsSurface<T>>(ref result2);
                return success;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                bool success = NurbsSurface3DHelper.RemoveKnot(
                    in surface,
                    removeKnot,
                    times,
                    direction,
                    out var result3
                );
                result = Unsafe.As<NurbsSurface<Vector3>, NurbsSurface<T>>(ref result3);
                return success;
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page209
        /// Algorithm A5.10
        /// Degree elevate a surface t times.
        /// </summary>
        /// <param name="times">Number of times to elevate the degree.</param>
        /// <param name="direction">The direction to elevate (UDirection or VDirection).</param>
        /// <returns>A new surface with elevated degree.</returns>
        public NurbsSurface<T> ElevateDegree(int times, SurfaceDirection direction)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.ElevateDegree(in surface, times, direction);
                return Unsafe.As<NurbsSurface<Vector2>, NurbsSurface<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.ElevateDegree(in surface, times, direction);
                return Unsafe.As<NurbsSurface<Vector3>, NurbsSurface<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page227
        /// Degree reduce U or V Direction Bezier-shape nurbs curve from degree to degree - 1.
        /// </summary>
        /// <param name="direction">The direction to reduce (UDirection or VDirection).</param>
        /// <param name="result">The resulting surface with reduced degree.</param>
        /// <returns>True if degree reduction was successful, false otherwise.</returns>
        public bool TryReduceDegree(SurfaceDirection direction, out NurbsSurface<T> result)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                bool success = NurbsSurface2DHelper.ReduceDegree(
                    in surface,
                    direction,
                    out var result2
                );
                result = Unsafe.As<NurbsSurface<Vector2>, NurbsSurface<T>>(ref result2);
                return success;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                bool success = NurbsSurface3DHelper.ReduceDegree(
                    in surface,
                    direction,
                    out var result3
                );
                result = Unsafe.As<NurbsSurface<Vector3>, NurbsSurface<T>>(ref result3);
                return success;
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page232
        /// Equally spaced parameter values on each candidate span.
        /// </summary>
        /// <param name="tessellatedPoints">Output: The tessellated surface points.</param>
        /// <param name="correspondingUVs">Output: The UV parameters corresponding to each point.</param>
        /// <param name="intervalsPerSpan">Number of intervals per knot span (default 100).</param>
        public void EquallyTessellate(
            out IReadOnlyList<T> tessellatedPoints,
            out IReadOnlyList<Vector2> correspondingUVs,
            int intervalsPerSpan = 100
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                NurbsSurface2DHelper.EquallyTessellate(
                    in surface,
                    out var points,
                    out correspondingUVs,
                    intervalsPerSpan
                );
                tessellatedPoints = Unsafe.As<IReadOnlyList<Vector2>, IReadOnlyList<T>>(ref points);
                return;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                NurbsSurface3DHelper.EquallyTessellate(
                    in surface,
                    out var points,
                    out correspondingUVs,
                    intervalsPerSpan
                );
                tessellatedPoints = Unsafe.As<IReadOnlyList<Vector3>, IReadOnlyList<T>>(ref points);
                return;
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Reparametrize the surface to a new parameter domain.
        /// The geometric shape remains unchanged, only the parameter values are rescaled.
        /// </summary>
        /// <param name="minU">New minimum U parameter.</param>
        /// <param name="maxU">New maximum U parameter.</param>
        /// <param name="minV">New minimum V parameter.</param>
        /// <param name="maxV">New maximum V parameter.</param>
        /// <returns>A new surface with rescaled parameter domain.</returns>
        public NurbsSurface<T> Reparametrize(
            float minU = 0f,
            float maxU = 1f,
            float minV = 0f,
            float maxV = 1f
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var result = NurbsSurface2DHelper.Reparametrize(in surface, minU, maxU, minV, maxV);
                return Unsafe.As<NurbsSurface<Vector2>, NurbsSurface<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var result = NurbsSurface3DHelper.Reparametrize(in surface, minU, maxU, minV, maxV);
                return Unsafe.As<NurbsSurface<Vector3>, NurbsSurface<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page232
        /// Point inversion:finding the corresponding parameter make S(u,v) = P.
        /// Find the UV parameter on the surface closest to the given point.
        /// Uses Newton-Raphson iteration with an initial guess from tessellation.
        /// </summary>
        /// <param name="givenPoint">The point to find the closest parameter for.</param>
        /// <param name="maxIterations">Maximum number of iterations (default 10).</param>
        /// <param name="tolerance">Convergence tolerance (default 1e-6).</param>
        /// <returns>The UV parameter closest to the given point.</returns>
        public Vector2 GetParamOnSurface(
            T givenPoint,
            int maxIterations = 10,
            double tolerance = 1e-6
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var point = Unsafe.As<T, Vector2>(ref givenPoint);
                return NurbsSurface2DHelper.GetParamOnSurface(
                    in surface,
                    point,
                    maxIterations,
                    tolerance
                );
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var point = Unsafe.As<T, Vector3>(ref givenPoint);
                return NurbsSurface3DHelper.GetParamOnSurface(
                    in surface,
                    point,
                    maxIterations,
                    tolerance
                );
            }

            throw new NotSupportedException(
                $"GetParamOnSurface is only supported for Vector3 surfaces."
            );
        }

        /// <summary>
        /// Find the UV parameter on the surface closest to the given point using
        /// a Geometric Surface Algorithm (GSA) that utilizes surface curvature.
        /// This method may converge better than Newton-Raphson for some cases.
        /// Experimental:
        /// According to https://jcst.ict.ac.cn/fileup/1000-9000/PDF/2019-6-9-9388.pdf
        /// A Geometric Strategy Algorithm for Orthogonal Projection onto a Parametric Surface
        /// </summary>
        /// <param name="givenPoint">The point to find the closest parameter for.</param>
        /// <param name="maxIterations">Maximum number of iterations (default 1000).</param>
        /// <param name="tolerance">Convergence tolerance (default 1e-10).</param>
        /// <returns>The UV parameter closest to the given point.</returns>
        public Vector2 GetParamOnSurfaceByGSA(
            T givenPoint,
            int maxIterations = 1000,
            double tolerance = 1e-10
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var point = Unsafe.As<T, Vector2>(ref givenPoint);
                return NurbsSurface2DHelper.GetParamOnSurfaceByGSA(
                    in surface,
                    point,
                    maxIterations,
                    tolerance
                );
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var point = Unsafe.As<T, Vector3>(ref givenPoint);
                return NurbsSurface3DHelper.GetParamOnSurfaceByGSA(
                    in surface,
                    point,
                    maxIterations,
                    tolerance
                );
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page235
        /// Surface Tangent Vector Inversion: finding the corresponding UV tangent [du dv] make T = Su*du+Sv*dv.
        /// Compute the UV tangent direction corresponding to a given tangent vector.
        /// Projects the tangent onto the surface's parameter space using the first fundamental form.
        /// </summary>
        /// <param name="param">The UV parameter on the surface.</param>
        /// <param name="tangent">The tangent vector in world space.</param>
        /// <param name="uvTangent">Output: The corresponding UV tangent direction.</param>
        /// <returns>True if the UV tangent was computed successfully, false if degenerate.</returns>
        public bool TryGetUVTangent(Vector2 param, T tangent, out Vector2 uvTangent)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                var tan = Unsafe.As<T, Vector2>(ref tangent);
                return NurbsSurface2DHelper.GetUVTangent(in surface, param, tan, out uvTangent);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                var tan = Unsafe.As<T, Vector3>(ref tangent);
                return NurbsSurface3DHelper.GetUVTangent(in surface, param, tan, out uvTangent);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Compute the approximate area of the surface using numerical integration.
        /// </summary>
        /// <param name="type">The integration method to use (Simpson, GaussLegendre, or Chebyshev).</param>
        /// <param name="tolerance">Tolerance for adaptive methods (default 1e-6).</param>
        /// <returns>The approximate surface area.</returns>
        public double ApproximateArea(
            IntegratorType type = IntegratorType.GaussLegendre,
            double tolerance = 1e-6
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                return NurbsSurface2DHelper.ApproximateArea(in surface, type, tolerance);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                return NurbsSurface3DHelper.ApproximateArea(in surface, type, tolerance);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <summary>
        /// Check if the surface is closed in the specified direction.
        /// A surface is closed in a direction if all iso-curves in that direction are closed.
        ///  [0][0]  [0][1] ... ...  [0][m]     ------- v direction
        ///  [1][0]  [1][1] ... ...  [1][m]    |
        ///    .                               |
        ///    .                               u direction
        ///    .
        ///  [n][0]  [n][1] ... ...  [n][m]
        /// </summary>
        /// <param name="direction">The direction to check (UDirection or VDirection).</param>
        /// <returns>True if the surface is closed in the specified direction.</returns>
        public bool IsClosed(SurfaceDirection direction)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector2>>(ref this);
                return NurbsSurface2DHelper.IsClosed(in surface, direction);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var surface = ref Unsafe.As<NurbsSurface<T>, NurbsSurface<Vector3>>(ref this);
                return NurbsSurface3DHelper.IsClosed(in surface, direction);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }
    }
}
