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
    public record struct BezierSurface<T>
        where T : struct
    {
        public int DegreeU { get; }
        public int DegreeV { get; }
        public IReadOnlyList<IReadOnlyList<ControlPoint<T>>> ControlPoints { get; }

        public BezierSurface(
            int degreeU,
            int degreeV,
            IReadOnlyList<IReadOnlyList<ControlPoint<T>>> controlPoints
        )
        {
            DegreeU = degreeU;
            DegreeV = degreeV;
            ControlPoints = controlPoints;

            Check();
        }

        public BezierSurface(
            int degreeU,
            int degreeV,
            IReadOnlyList<IReadOnlyList<T>> controlPoints
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

            Check();
        }

        public BezierSurface(Int2 controlPointsCount, IReadOnlyList<IReadOnlyList<T>> controlPoints)
        {
            DegreeU = controlPointsCount.X - 1;
            DegreeV = controlPointsCount.Y - 1;

            ControlPoints = controlPoints
                .Select(row =>
                    (IReadOnlyList<ControlPoint<T>>)
                        row.Select(cp => new ControlPoint<T>(cp)).ToImmutableArray()
                )
                .ToImmutableArray();

            Check();
        }

        public void Check()
        {
            var degreeU = DegreeU;
            var degreeV = DegreeV;

            Validate.Argument(degreeU > 0, nameof(DegreeU), "DegreeU must be greater than zero.");
            Validate.Argument(degreeV > 0, nameof(DegreeV), "DegreeV must be greater than zero.");
            Validate.Argument(
                ControlPoints.Count > 0,
                nameof(ControlPoints),
                "ControlPoints must contain one point at least."
            );
            Validate.Argument(
                ControlPoints.Count == degreeU + 1,
                nameof(ControlPoints),
                "ControlPoints row count must be equal to DegreeU + 1."
            );
            Validate.Argument(
                ControlPoints.All(row => row.Count == degreeV + 1),
                nameof(ControlPoints),
                "ControlPoints column count must be equal to DegreeV + 1."
            );
        }

        /// <summary>
        /// The NURBS Book 2nd Edition Page39
        /// Algorithm A1.7
        /// Compute a point on a Bezier surface by the deCasteljau.
        ///
        /// Controlpoints with (n+1) rows * (m+1) columns
        ///
        ///  [0][0]  [0][1] ... ...  [0][m]     ------- v direction
        ///  [1][0]  [1][1] ... ...  [1][m]    |
        ///    .                               |
        ///    .                               u direction
        ///    .
        ///  [n][0]  [n][1] ... ...  [n][m]
        ///
        /// Rational Bezier Surface:Use XYZW
        /// </summary>
        public T GetPointOnSurfaceByDeCasteljau(Vector2 uv)
        {
            Validate.Argument(
                uv.X >= 0 && uv.X <= 1,
                nameof(uv),
                "The u parameter must be in the range [0,1]."
            );

            Validate.Argument(
                uv.Y >= 0 && uv.Y <= 1,
                nameof(uv),
                "The v parameter must be in the range [0,1]."
            );

            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<BezierSurface<T>, BezierSurface<Vector2>>(ref this);
                var result = BezierSurface2DHelper.GetPointOnSurfaceByDeCasteljau(in curve2, uv);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<BezierSurface<T>, BezierSurface<Vector3>>(ref this);
                var result = BezierSurface3DHelper.GetPointOnSurfaceByDeCasteljau(in curve3, uv);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }

        /// <inheritdoc cref="GetPointOnSurfaceByDeCasteljau(Vector2)"/>
        public T GetPointOnSurface(Vector2 uv)
        {
            return GetPointOnSurfaceByDeCasteljau(uv);
        }
    }
}
