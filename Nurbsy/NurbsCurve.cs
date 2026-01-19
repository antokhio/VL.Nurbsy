/*
 * Original Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * Ported by:
 * Anton Kalabukhov (antokhio)
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using System.Collections.Immutable;
using System.Runtime.CompilerServices;
using Nurbsy.Algorithm;
using Nurbsy.Helpers;
using Stride.Core.Mathematics;

namespace Nurbsy
{
    /// <summary>
    /// Adaptive base class for NURBS curves.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public record struct NurbsCurve<T>
        where T : struct
    {
        /// <summary>
        /// Gets the degree of the NURBS curve.
        /// </summary>
        public int Degree { get; }

        /// <summary>
        /// Gets the control points of the NURBS curve.
        /// </summary>
        public IReadOnlyList<ControlPoint<T>> ControlPoints { get; }

        /// <summary>
        /// Gets the knot vector of the NURBS curve.
        /// </summary>
        public IReadOnlyList<double> Knots { get; }

        /// <summary>
        /// Initializes a new instance of the <see cref="NurbsCurve{T}"/> struct with clamped knots generated automatically.
        /// </summary>
        /// <param name="degree">The degree of the curve. Must be greater than zero.</param>
        /// <param name="controlPoints">The control points of the curve.</param>
        public NurbsCurve(int degree, IReadOnlyList<ControlPoint<T>> controlPoints)
        {
            Degree = degree;
            ControlPoints = controlPoints;
            Knots = KnotsUtils.GenerateClampedKnots(Degree, ControlPoints.Count);

            Check();
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="NurbsCurve{T}"/> struct with a specified knot vector.
        /// </summary>
        /// <param name="degree">The degree of the curve. Must be greater than zero.</param>
        /// <param name="controlPoints">The control points of the curve.</param>
        /// <param name="knots">The knot vector. Must be non-decreasing.</param>
        public NurbsCurve(
            int degree,
            IReadOnlyList<ControlPoint<T>> controlPoints,
            IReadOnlyList<double> knots
        )
        {
            Degree = degree;
            ControlPoints = controlPoints;
            Knots = knots;

            Check();
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="NurbsCurve{T}"/> struct with raw points as control points (weight 1.0).
        /// </summary>
        /// <param name="degree">The degree of the curve. Must be greater than zero.</param>
        /// <param name="controlPoints">The points to be used as control points.</param>
        /// <param name="knots">The knot vector.</param>
        public NurbsCurve(int degree, IReadOnlyList<T> controlPoints, IReadOnlyList<double> knots)
        {
            Degree = degree;
            ControlPoints = controlPoints.Select(cp => new ControlPoint<T>(cp)).ToImmutableArray();
            Knots = knots;

            Check();
        }

        /// <summary>
        /// Check curve whether fits NURBS.
        /// </summary>
        /// <exception cref="ArgumentException">Thrown if degree is non-positive, knots/control points are empty, knots are invalid, or dimensions do not match (m = n + p + 1).</exception>
        public void Check()
        {
            Validate.Argument(Degree > 0, nameof(Degree), "Degree must be gretater then zero.");
            Validate.Argument(
                Knots.Count > 0,
                nameof(Knots),
                "Knot vector must contain at least one knot."
            );
            Validate.Argument(
                Validate.IsValidKnots(Knots),
                nameof(Knots),
                "Knot must be a nondecreasing sequence of real numbers."
            );
            Validate.Argument(
                ControlPoints.Count > 0,
                nameof(ControlPoints),
                "ControlPoints must contain one point at least."
            );
            Validate.Argument(
                Validate.IsValidNURBS(Degree, ControlPoints.Count, Knots.Count),
                nameof(NurbsCurve<T>),
                "Arguments must be fit: m = n + p + 1"
            );
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 124, Algorithm A4.1.
        /// Compute point on rational B-spline curve.
        /// </summary>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>The point on the curve.</returns>
        public T GetPointOnCurve(float paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.GetPointOnCurve(curve2, paramT);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.GetPointOnCurve(curve3, paramT);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 155, Algorithm A5.2.
        /// Computes point on rational B-spline curve using the corner cutting algorithm.
        /// </summary>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>The point on the curve.</returns>
        public T GetPointOnCurveByCornerCut(float paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.GetPointOnCurveByCornerCut(curve2, paramT);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.GetPointOnCurveByCornerCut(curve3, paramT);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 127, Algorithm A4.2.
        /// Compute C(paramT) derivatives from Cw(paramT) derivatives.
        /// </summary>
        /// <param name="derivative">The number of derivatives to compute.</param>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>A list of derivatives where index 0 is the point itself.</returns>
        public IReadOnlyList<T> GetDerivatives(int derivative, float paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.ComputeRationalCurveDerivatives(
                    curve2,
                    derivative,
                    paramT
                );
                return Unsafe.As<IReadOnlyList<Vector2>, IReadOnlyList<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.ComputeRationalCurveDerivatives(
                    curve3,
                    derivative,
                    paramT
                );
                return Unsafe.As<IReadOnlyList<Vector3>, IReadOnlyList<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Computer left and right hand derivatives.
        /// Checks if derivatives can be computed at the specified parameter (considering knot multiplicity).
        /// </summary>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>True if derivatives can be computed; otherwise, false.</returns>
        public bool CanComputeDerivative(float paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.CanComputeDerivative(curve2, paramT);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.CanComputeDerivative(curve3, paramT);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Calculate curve curvature.
        /// </summary>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>The curvature value.</returns>
        public float GetCurvature(float paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return (float)NurbsCurve2DHelper.GetCurvature(curve2, paramT);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return (float)NurbsCurve3DHelper.GetCurvature(curve3, paramT);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Calculate curve torsion.
        /// </summary>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>The torsion value.</returns>
        public float GetTorsion(float paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return (float)NurbsCurve2DHelper.GetTorsion(curve2, paramT);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return (float)NurbsCurve3DHelper.GetTorsion(curve3, paramT);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 479.
        /// Calculate curve normal direction.
        /// </summary>
        /// <param name="normalType">The type of normal to calculate (e.g., Normal, Binormal).</param>
        /// <param name="paramT">The parameter value.</param>
        /// <returns>The normal vector.</returns>
        public T GetNormal(CurveNormal normalType, float paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.GetNormal(curve2, normalType, paramT);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.GetNormal(curve3, normalType, paramT);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 481.
        /// Projection normal method invented by Siltanen and Woodward.
        /// </summary>
        /// <returns>A list of projected normal vectors.</returns>
        public IReadOnlyList<T> ProjectNormal()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.ProjectNormal(curve2);
                return Unsafe.As<IReadOnlyList<Vector2>, IReadOnlyList<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.ProjectNormal(curve3);
                return Unsafe.As<IReadOnlyList<Vector3>, IReadOnlyList<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 173, Algorithm A5.6.
        /// Decompose curve into Bezier segements.
        /// </summary>
        /// <returns>A list of BezierCurve{T}.</returns>
        public IReadOnlyList<BezierCurve<T>> DecomposeToBeziers()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.DecomposeToBeziers(curve2);
                return Unsafe.As<
                    IReadOnlyList<BezierCurve<Vector2>>,
                    IReadOnlyList<BezierCurve<T>>
                >(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.DecomposeToBeziers(curve3);
                return Unsafe.As<
                    IReadOnlyList<BezierCurve<Vector3>>,
                    IReadOnlyList<BezierCurve<T>>
                >(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 230.
        /// Equally spaced parameter values on each candidate span.
        /// </summary>
        /// <returns>A tuple containing the tessellated points and their corresponding knot parameters.</returns>
        public (
            IReadOnlyList<T> TessellatedPoints,
            IReadOnlyList<double> CorrespondingKnots
        ) EquallyTessellate()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var (pts, knots) = NurbsCurve2DHelper.EquallyTessellate(curve2);
                var ptsT = Unsafe.As<IReadOnlyList<Vector2>, IReadOnlyList<T>>(ref pts);
                return (ptsT, knots);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var (pts, knots) = NurbsCurve3DHelper.EquallyTessellate(curve3);
                var ptsT = Unsafe.As<IReadOnlyList<Vector3>, IReadOnlyList<T>>(ref pts);
                return (ptsT, knots);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Tessellate nurbs curve.
        /// </summary>
        /// <param name="tolerance">The tolerance for tessellation.</param>
        /// <returns>A list of tessellated points.</returns>
        public IReadOnlyList<T> Tessellate(float tolerance = 0.0001f)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.Tessellate(curve2, tolerance);
                return Unsafe.As<IReadOnlyList<Vector2>, IReadOnlyList<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.Tessellate(curve3, tolerance);
                return Unsafe.As<IReadOnlyList<Vector3>, IReadOnlyList<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 230.
        /// Point inversion: finding the corresponding parameter make C(u) = P.
        /// </summary>
        /// <param name="point">The point on or near the curve.</param>
        /// <returns>The parameter value on the curve.</returns>
        public double GetParamOnCurve(T point)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                ref var p = ref Unsafe.As<T, Vector2>(ref point);
                return NurbsCurve2DHelper.GetParamOnCurve(curve2, p);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                ref var p = ref Unsafe.As<T, Vector3>(ref point);
                return NurbsCurve3DHelper.GetParamOnCurve(curve3, p);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 255.
        /// Reparameterization using a linear rational function : (alpha * u + beta)/(gamma * u + delta).
        /// </summary>
        /// <param name="alpha">Coefficient alpha.</param>
        /// <param name="beta">Coefficient beta.</param>
        /// <param name="gamma">Coefficient gamma.</param>
        /// <param name="delta">Coefficient delta.</param>
        /// <returns>A reparameterized NURBS curve.</returns>
        public NurbsCurve<T> Reparametrize(double alpha, double beta, double gamma, double delta)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.Reparametrize(curve2, alpha, beta, gamma, delta);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.Reparametrize(curve3, alpha, beta, gamma, delta);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 241.
        /// Reparameterization of curve values linearly to a new range.
        /// </summary>
        /// <param name="min">The new start parameter.</param>
        /// <param name="max">The new end parameter.</param>
        /// <returns>A reparameterized NURBS curve.</returns>
        public NurbsCurve<T> Reparametrize(double min, double max)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.Reparametrize(curve2, min, max);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.Reparametrize(curve3, min, max);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 164, Algorithm A5.4.
        /// Refine curve knot vector.
        /// </summary>
        /// <param name="insertKnotElements">The knots to insert.</param>
        /// <returns>A new NURBS curve with the refined knot vector.</returns>
        public NurbsCurve<T> RefineKnotVector(IReadOnlyList<double> insertKnotElements)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var refined2 = NurbsCurve2DHelper.RefineKnotVector(curve2, insertKnotElements);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref refined2);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var refined3 = NurbsCurve3DHelper.RefineKnotVector(curve3, insertKnotElements);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref refined3);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 263.
        /// Curve reverse, but not use reparameterization.
        /// </summary>
        /// <returns>The reversed NURBS curve.</returns>
        public NurbsCurve<T> Reverse()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.Reverse(curve2);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.Reverse(curve3);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Split curve at certain parameter.
        /// </summary>
        /// <param name="parameter">The parameter at which to split the curve.</param>
        /// <param name="left">The resulting left part of the curve.</param>
        /// <param name="right">The resulting right part of the curve.</param>
        /// <returns>True if the split was successful; otherwise, false.</returns>
        public bool SplitAt(double parameter, out NurbsCurve<T> left, out NurbsCurve<T> right)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.SplitAt(curve2, parameter, out var l2, out var r2))
                {
                    left = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref l2);
                    right = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref r2);
                    return true;
                }
                left = default;
                right = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.SplitAt(curve3, parameter, out var l3, out var r3))
                {
                    left = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref l3);
                    right = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref r3);
                    return true;
                }
                left = default;
                right = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Segment curve.
        /// </summary>
        /// <param name="startParameter">The start parameter of the segment.</param>
        /// <param name="endParameter">The end parameter of the segment.</param>
        /// <param name="segment">The extracted curve segment.</param>
        /// <returns>True if the segmentation was successful; otherwise, false.</returns>
        public bool Segment(double startParameter, double endParameter, out NurbsCurve<T> segment)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.Segment(curve2, startParameter, endParameter, out var res2))
                {
                    segment = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                segment = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.Segment(curve3, startParameter, endParameter, out var res3))
                {
                    segment = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                segment = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 151, Algorithm A5.1.
        /// Curve knot insertion.
        /// Note that multiplicity + times must be &lt;= degree.
        /// </summary>
        /// <param name="insertKnot">The knot value to insert.</param>
        /// <param name="times">The number of times to insert the knot.</param>
        /// <param name="result">The resulting curve with inserted knots.</param>
        /// <returns>The number of knots successfully inserted.</returns>
        public int InsertKnot(double insertKnot, int times, out NurbsCurve<T> result)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                int inserted = NurbsCurve2DHelper.InsertKnot(
                    curve2,
                    insertKnot,
                    times,
                    out var res2
                );
                result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                return inserted;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                int inserted = NurbsCurve3DHelper.InsertKnot(
                    curve3,
                    insertKnot,
                    times,
                    out var res3
                );
                result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                return inserted;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 185, Algorithm A5.8.
        /// Curve knot removal.
        /// </summary>
        /// <param name="removeKnot">The knot value to remove.</param>
        /// <param name="times">The number of times to remove the knot.</param>
        /// <param name="result">The resulting curve with removed knots.</param>
        /// <returns>True if the removal was successful; otherwise, false.</returns>
        public bool RemoveKnot(double removeKnot, int times, out NurbsCurve<T> result)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.RemoveKnot(curve2, removeKnot, times, out var res2))
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default; // Should be curve copy?
                return false;
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.RemoveKnot(curve3, removeKnot, times, out var res3))
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Remove excessive knots.
        /// </summary>
        /// <returns>The curve without excessive knots.</returns>
        public NurbsCurve<T> RemoveExcessiveKnots()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.RemoveExcessiveKnots(curve2);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.RemoveExcessiveKnots(curve3);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 206, Algorithm A5.9.
        /// Degree elevate a curve t times.
        /// </summary>
        /// <param name="times">Number of times to elevate degree.</param>
        /// <returns>Degree elevated curve.</returns>
        public NurbsCurve<T> ElevateDegree(int times)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.ElevateDegree(curve2, times);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.ElevateDegree(curve3, times);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 223, Algorithm A5.11.
        /// Degree reduce a bezier-shape nurbs curve from degree to degree - 1.
        /// </summary>
        /// <param name="result">The resulting curve with reduced degree.</param>
        /// <returns>True if reduction successful, otherwise false.</returns>
        public bool ReduceDegree(out NurbsCurve<T> result)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.ReduceDegree(curve2, out var res2))
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default;
                return false;
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.ReduceDegree(curve3, out var res3))
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 429, Algorithm A9.9.
        /// Remove knots from curve by given bound.
        /// </summary>
        /// <param name="parameters">Parameters to consider.</param>
        /// <param name="errors">List to store errors.</param>
        /// <param name="maxError">Maximum error allowed.</param>
        /// <param name="result">The resulting curve.</param>
        /// <returns>True if knots removed within bound, otherwise false.</returns>
        public bool RemoveKnotsByGivenBound(
            IReadOnlyList<double> parameters,
            List<double> errors,
            double maxError,
            out NurbsCurve<T> result
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (
                    NurbsCurve2DHelper.RemoveKnotsByGivenBound(
                        curve2,
                        parameters,
                        errors,
                        maxError,
                        out var res2
                    )
                )
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (
                    NurbsCurve3DHelper.RemoveKnotsByGivenBound(
                        curve3,
                        parameters,
                        errors,
                        maxError,
                        out var res3
                    )
                )
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Merge two connected curves to one curve.
        /// </summary>
        /// <param name="left">The left curve.</param>
        /// <param name="right">The right curve.</param>
        /// <param name="result">The resulting merged curve.</param>
        /// <returns>True if merge successful, otherwise false.</returns>
        public static bool Merge(NurbsCurve<T> left, NurbsCurve<T> right, out NurbsCurve<T> result)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var l2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref left);
                ref var r2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref right);
                if (NurbsCurve2DHelper.Merge(l2, r2, out var res2))
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var l3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref left);
                ref var r3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref right);
                if (NurbsCurve3DHelper.Merge(l3, r3, out var res3))
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 511.
        /// Reposition an arbitrary control point.
        /// </summary>
        /// <param name="parameter">The parameter on curve.</param>
        /// <param name="moveIndex">Index of control point to move.</param>
        /// <param name="moveDirection">Direction to move.</param>
        /// <param name="moveDistance">Distance to move.</param>
        /// <param name="result">The resulting curve.</param>
        /// <returns>True if reposition successful, otherwise false.</returns>
        public bool ControlPointReposition(
            double parameter,
            int moveIndex,
            T moveDirection,
            double moveDistance,
            out NurbsCurve<T> result
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                ref var dir2 = ref Unsafe.As<T, Vector2>(ref moveDirection);
                if (
                    NurbsCurve2DHelper.ControlPointReposition(
                        curve2,
                        parameter,
                        moveIndex,
                        dir2,
                        moveDistance,
                        out var res2
                    )
                )
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                ref var dir3 = ref Unsafe.As<T, Vector3>(ref moveDirection);
                if (
                    NurbsCurve3DHelper.ControlPointReposition(
                        curve3,
                        parameter,
                        moveIndex,
                        dir3,
                        moveDistance,
                        out var res3
                    )
                )
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 526.
        /// Modify two neighboring curve weights. (moveIndex and moveIndex + 1).
        /// </summary>
        /// <param name="parameter">Parameter on curve.</param>
        /// <param name="moveIndex">Index of first weight.</param>
        /// <param name="moveDistance">Distance to move.</param>
        /// <param name="scale">Scale factor.</param>
        /// <param name="result">Resulting curve.</param>
        /// <returns>True if modification successful, otherwise false.</returns>
        public bool NeighborWeightsModification(
            double parameter,
            int moveIndex,
            double moveDistance,
            double scale,
            out NurbsCurve<T> result
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (
                    NurbsCurve2DHelper.NeighborWeightsModification(
                        curve2,
                        parameter,
                        moveIndex,
                        moveDistance,
                        scale,
                        out var res2
                    )
                )
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (
                    NurbsCurve3DHelper.NeighborWeightsModification(
                        curve3,
                        parameter,
                        moveIndex,
                        moveDistance,
                        scale,
                        out var res3
                    )
                )
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 533.
        /// Warping the curve.
        /// </summary>
        /// <param name="warpShape">Shape of warp.</param>
        /// <param name="warpDistance">Distance of warp.</param>
        /// <param name="planeNormal">Normal of plane.</param>
        /// <param name="startParameter">Start parameter of warp.</param>
        /// <param name="endParameter">End parameter of warp.</param>
        /// <returns>Warped curve.</returns>
        public NurbsCurve<T> Warping(
            IReadOnlyList<double> warpShape,
            double warpDistance,
            T planeNormal,
            double startParameter,
            double endParameter
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.Warping(
                    curve2,
                    warpShape,
                    warpDistance,
                    startParameter,
                    endParameter
                );
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                ref var pn3 = ref Unsafe.As<T, Vector3>(ref planeNormal);
                var result = NurbsCurve3DHelper.Warping(
                    curve3,
                    warpShape,
                    warpDistance,
                    pn3,
                    startParameter,
                    endParameter
                );
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 542.
        /// Flattening the curve.
        /// </summary>
        /// <param name="lineStart">Line start point.</param>
        /// <param name="lineEnd">Line end point.</param>
        /// <param name="startParameter">Start parameter.</param>
        /// <param name="endParameter">End parameter.</param>
        /// <param name="result">Flattened curve.</param>
        /// <returns>True if flattening successful, otherwise false.</returns>
        public bool Flattening(
            T lineStart,
            T lineEnd,
            double startParameter,
            double endParameter,
            out NurbsCurve<T> result
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                ref var start2 = ref Unsafe.As<T, Vector2>(ref lineStart);
                ref var end2 = ref Unsafe.As<T, Vector2>(ref lineEnd);
                if (
                    NurbsCurve2DHelper.Flattening(
                        curve2,
                        start2,
                        end2,
                        startParameter,
                        endParameter,
                        out var res2
                    )
                )
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                ref var start3 = ref Unsafe.As<T, Vector3>(ref lineStart);
                ref var end3 = ref Unsafe.As<T, Vector3>(ref lineEnd);
                if (
                    NurbsCurve3DHelper.Flattening(
                        curve3,
                        start3,
                        end3,
                        startParameter,
                        endParameter,
                        out var res3
                    )
                )
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 547.
        /// Bending the curve.
        /// </summary>
        /// <param name="startParameter">Start parameter.</param>
        /// <param name="endParameter">End parameter.</param>
        /// <param name="bendCenter">Center of bend.</param>
        /// <param name="radius">Radius of bend.</param>
        /// <param name="crossRatio">Cross ratio.</param>
        /// <returns>Bended curve.</returns>
        public NurbsCurve<T> Bending(
            double startParameter,
            double endParameter,
            T bendCenter,
            double radius,
            double crossRatio
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                ref var center2 = ref Unsafe.As<T, Vector2>(ref bendCenter);
                var result = NurbsCurve2DHelper.Bending(
                    curve2,
                    startParameter,
                    endParameter,
                    center2,
                    radius,
                    crossRatio
                );
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                ref var center3 = ref Unsafe.As<T, Vector3>(ref bendCenter);
                var result = NurbsCurve3DHelper.Bending(
                    curve3,
                    startParameter,
                    endParameter,
                    center3,
                    radius,
                    crossRatio
                );
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 572.
        /// Clamp a unclamped curve.
        /// </summary>
        /// <returns>Clamped curve.</returns>
        public NurbsCurve<T> ToClampCurve()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.ToClampCurve(curve2);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.ToClampCurve(curve3);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 577, Algorithm A12.1.
        /// Unclamp a clamped curve.
        /// </summary>
        /// <returns>Unclamped curve.</returns>
        public NurbsCurve<T> ToUnclampCurve()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.ToUnclampCurve(curve2);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.ToUnclampCurve(curve3);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Calculate curve arc length.
        /// </summary>
        /// <param name="type">Integrator type to use.</param>
        /// <returns>Approximate length of curve.</returns>
        public double ApproximateLength(IntegratorType type)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.ApproximateLength(curve2, type);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.ApproximateLength(curve3, type);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Calculate parameter makes first segment length equals to given length.
        /// </summary>
        /// <param name="lenght">Given Lenght</param>
        /// <param name="type">IntegratorType</param>
        /// <returns>Parameter on curve.</returns>
        /// <exception cref="NotSupportedException"></exception>
        public double GetParamOnCurve(float lenght, IntegratorType type)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.GetParamOnCurve(curve2, lenght, type);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.GetParamOnCurve(curve3, lenght, type);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Calculate parameters makes every segments length equals to given length.
        /// </summary>
        /// <param name="length">Given Lenght</param>
        /// <param name="type">IntegratorType</param>
        /// <returns>IReadOnlyList Vector2 or Vecotor3</returns>
        /// <exception cref="NotSupportedException">Adaptive for Vector2 or Vector3</exception>
        public IReadOnlyList<double> GetParamsOnCurve(
            float length = 0.0001f,
            IntegratorType type = IntegratorType.GaussLegendre
        )
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.GetParamsOnCurve(curve2, length, type);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.GetParamsOnCurve(curve3, length, type);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 575.
        /// UnClamped, uniform and C[degree-1] continues closed curve.
        /// </summary>
        /// <returns>True if periodic, otherwise false.</returns>
        public bool IsPeriodic()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.IsPeriodic(curve2);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.IsPeriodic(curve3);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page 572.
        /// Detemine curve is clamped.
        /// </summary>
        /// <returns>True if clamped, otherwise false.</returns>
        public bool IsClamp()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.IsClamp(curve2);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.IsClamp(curve3);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Detemine curve whether is linear.
        /// </summary>
        /// <returns>True if linear, otherwise false.</returns>
        public bool IsLinear()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.IsLinear(curve2);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.IsLinear(curve3);
            }

            return false;
        }

        /// <summary>
        /// Detemine curve whether is arc.
        /// </summary>
        /// <param name="center">Center of arc.</param>
        /// <param name="radius">Radius of arc.</param>
        /// <returns>True if arc, otherwise false.</returns>
        public bool IsArc(out T center, out double radius)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.IsArc(curve2, out var c2, out radius))
                {
                    center = Unsafe.As<Vector2, T>(ref c2);
                    return true;
                }
                center = default;
                radius = 0.0;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.IsArc(curve3, out var c3, out radius))
                {
                    center = Unsafe.As<Vector3, T>(ref c3);
                    return true;
                }
                center = default;
                radius = 0.0;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        /// <summary>
        /// Detemine curve is closed.
        /// Close means end point equals start point or points overlap.
        /// </summary>
        /// <returns>True if closed, otherwise false.</returns>
        public bool IsClosed()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.IsClosed(curve2);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.IsClosed(curve3);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
