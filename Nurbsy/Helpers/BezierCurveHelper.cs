/*
 * Ported by:
 * 2025 - Anton Kalabuhov (antokhio)
 *
 * Original Author:
 * 2023 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 * Based on LNLib: https://github.com/BIMCoderLiang/LNLib
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using System.Runtime.CompilerServices;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    public static class BezierCurveHelper
    {
        public static T GetPointOnCurveByDeCasteljau<T>(BezierCurve<T> curve, double paramT)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<BezierCurve<T>, BezierCurve<Vector2>>(ref curve);
                var result = BezierCurve2DHelper.GetPointOnCurveByDeCasteljau(curve2, paramT);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<BezierCurve<T>, BezierCurve<Vector3>>(ref curve);
                var result = BezierCurve3DHelper.GetPointOnCurveByDeCasteljau(curve3, paramT);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static T GetPointOnCurveByBernstein<T>(BezierCurve<T> curve, double paramT)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<BezierCurve<T>, BezierCurve<Vector2>>(ref curve);
                var result = BezierCurve2DHelper.GetPointOnCurveByBernstein(curve2, paramT);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<BezierCurve<T>, BezierCurve<Vector3>>(ref curve);
                var result = BezierCurve3DHelper.GetPointOnCurveByBernstein(curve3, paramT);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static T GetPointOnQuadraticArc<T>(
            ControlPoint<T> startPoint,
            ControlPoint<T> middlePoint,
            ControlPoint<T> endPoint,
            double paramT
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var p0 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector2>>(ref startPoint);
                ref var p1 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector2>>(ref middlePoint);
                ref var p2 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector2>>(ref endPoint);

                var result = BezierCurve2DHelper.GetPointOnQuadraticArc(
                    ref p0,
                    ref p1,
                    ref p2,
                    paramT
                );
                return Unsafe.As<Vector2, T>(ref result);
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var p0 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector3>>(ref startPoint);
                ref var p1 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector3>>(ref middlePoint);
                ref var p2 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector3>>(ref endPoint);

                var result = BezierCurve3DHelper.GetPointOnQuadraticArc(
                    ref p0,
                    ref p1,
                    ref p2,
                    paramT
                );
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static IReadOnlyList<ControlPoint<T>>? ComputeMiddleControlPointsOnQuadraticCurve<T>(
            T startPoint,
            T startTangent,
            T endPoint,
            T endTangent
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var sP = ref Unsafe.As<T, Vector2>(ref startPoint);
                ref var sT = ref Unsafe.As<T, Vector2>(ref startTangent);
                ref var eP = ref Unsafe.As<T, Vector2>(ref endPoint);
                ref var eT = ref Unsafe.As<T, Vector2>(ref endTangent);

                IList<ControlPoint<Vector2>> list = new List<ControlPoint<Vector2>>();
                if (
                    BezierCurve2DHelper.ComputeMiddleControlPointsOnQuadraticCurve(
                        ref sP,
                        ref sT,
                        ref eP,
                        ref eT,
                        ref list
                    )
                )
                {
                    var array = list.ToArray();
                    return Unsafe.As<ControlPoint<Vector2>[], ControlPoint<T>[]>(ref array);
                }
                return null;
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var sP = ref Unsafe.As<T, Vector3>(ref startPoint);
                ref var sT = ref Unsafe.As<T, Vector3>(ref startTangent);
                ref var eP = ref Unsafe.As<T, Vector3>(ref endPoint);
                ref var eT = ref Unsafe.As<T, Vector3>(ref endTangent);

                IList<ControlPoint<Vector3>> list = new List<ControlPoint<Vector3>>();
                if (
                    BezierCurve3DHelper.ComputeMiddleControlPointsOnQuadraticCurve(
                        ref sP,
                        ref sT,
                        ref eP,
                        ref eT,
                        ref list
                    )
                )
                {
                    var array = list.ToArray();
                    return Unsafe.As<ControlPoint<Vector3>[], ControlPoint<T>[]>(ref array);
                }
                return null;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
