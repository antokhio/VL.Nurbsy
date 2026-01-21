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

namespace Nurbsy.Algorithm
{
    public static class Intersection
    {
        public static CurveCurveIntersectionType ComputeRays<T>(
            T point0,
            T vector0,
            T point1,
            T vector1,
            out double param0,
            out double param1,
            out T intersectPoint
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var p0 = Unsafe.As<T, Vector2>(ref point0);
                var v0 = Unsafe.As<T, Vector2>(ref vector0);
                var p1 = Unsafe.As<T, Vector2>(ref point1);
                var v1 = Unsafe.As<T, Vector2>(ref vector1);

                var result = Intersection2D.ComputeRays(
                    p0,
                    v0,
                    p1,
                    v1,
                    out param0,
                    out param1,
                    out var intersect
                );
                intersectPoint = Unsafe.As<Vector2, T>(ref intersect);
                return result;
            }
            else if (typeof(T) == typeof(Vector3))
            {
                var p0 = Unsafe.As<T, Vector3>(ref point0);
                var v0 = Unsafe.As<T, Vector3>(ref vector0);
                var p1 = Unsafe.As<T, Vector3>(ref point1);
                var v1 = Unsafe.As<T, Vector3>(ref vector1);

                var result = Intersection3D.ComputeRays(
                    p0,
                    v0,
                    p1,
                    v1,
                    out param0,
                    out param1,
                    out var intersect
                );
                intersectPoint = Unsafe.As<Vector3, T>(ref intersect);
                return result;
            }

            throw new NotSupportedException(
                $"Type {typeof(T).Name} is not supported for ComputeRays."
            );
        }

        public static LinePlaneIntersectionType ComputeLineAndPlane<T>(
            T normal,
            T pointOnPlane,
            T pointOnLine,
            T lineDirection,
            out T intersectPoint
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector3))
            {
                var n = Unsafe.As<T, Vector3>(ref normal);
                var pp = Unsafe.As<T, Vector3>(ref pointOnPlane);
                var pl = Unsafe.As<T, Vector3>(ref pointOnLine);
                var ld = Unsafe.As<T, Vector3>(ref lineDirection);

                var result = Intersection3D.ComputeLineAndPlane(n, pp, pl, ld, out var intersect);
                intersectPoint = Unsafe.As<Vector3, T>(ref intersect);
                return result;
            }

            throw new NotSupportedException(
                $"Type {typeof(T).Name} is not supported for ComputeLineAndPlane."
            );
        }
    }
}
