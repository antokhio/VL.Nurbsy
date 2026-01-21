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
    internal struct Constants
    {
        public const double DoubleEpsilon = 1E-6;
        public const double DistanceEpsilon = 1E-4;
        public const double AngleEpsilon = 1E-2;
        public const double MaxDistance = 1E9;
        public const double Pi = 3.14159265358979323846;
        public const int NURBSMaxDegree = 7;
    }
}
