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

using System.Numerics;
using System.Runtime.CompilerServices;

namespace Nurbsy.Helpers
{
    internal static class ControlPointsHelper
    {
        public static ControlPoint<T> BlendControlPoints<T>(
            ControlPoint<T> cp0,
            ControlPoint<T> cp1,
            double alpha
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var c0 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector2>>(
                    ref Unsafe.AsRef(in cp0)
                );
                ref var c1 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector2>>(
                    ref Unsafe.AsRef(in cp1)
                );

                double w0 = c0.Weight;
                double w1 = c1.Weight;

                double newWeight = (1.0 - alpha) * w0 + alpha * w1;

                if (Math.Abs(newWeight) < 1e-9) // Handle singularity
                {
                    var zero = new ControlPoint<Vector2>(Vector2.Zero, 0.0);
                    return Unsafe.As<ControlPoint<Vector2>, ControlPoint<T>>(ref zero);
                }

                // Weighted interpolation: (w0 * v0 * (1-a) + w1 * v1 * a) / newWeight
                float k0 = (float)((1.0 - alpha) * w0);
                float k1 = (float)(alpha * w1);

                Vector2 newValue = (k0 * c0.Value + k1 * c1.Value) / (float)newWeight;

                // 3. Construct result and cast back to generic T
                var result = new ControlPoint<Vector2>(newValue, newWeight);
                return Unsafe.As<ControlPoint<Vector2>, ControlPoint<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var c0 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector3>>(
                    ref Unsafe.AsRef(in cp0)
                );
                ref var c1 = ref Unsafe.As<ControlPoint<T>, ControlPoint<Vector3>>(
                    ref Unsafe.AsRef(in cp1)
                );

                double w0 = c0.Weight;
                double w1 = c1.Weight;
                double newWeight = (1.0 - alpha) * w0 + alpha * w1;

                if (Math.Abs(newWeight) < 1e-9)
                {
                    var zero = new ControlPoint<Vector3>(Vector3.Zero, 0.0);
                    return Unsafe.As<ControlPoint<Vector3>, ControlPoint<T>>(ref zero);
                }

                float k0 = (float)((1.0 - alpha) * w0);
                float k1 = (float)(alpha * w1);

                Vector3 newValue = (k0 * c0.Value + k1 * c1.Value) / (float)newWeight;

                var result = new ControlPoint<Vector3>(newValue, newWeight);
                return Unsafe.As<ControlPoint<Vector3>, ControlPoint<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T)} is not supported.");
        }
    }
}
