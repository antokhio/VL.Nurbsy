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
namespace Nurbsy
{
    public record struct ControlPoint<T>
        where T : struct
    {
        public T Value { get; set; }
        public double Weight { get; set; }

        public ControlPoint(T value)
        {
            Value = value;
            Weight = 1f;
        }

        public ControlPoint(T value, double weight)
        {
            Value = value;
            Weight = weight;
        }
    }
}
