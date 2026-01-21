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

using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    public static class BezierSurface2DHelper
    {
        public static Vector2 GetPointOnSurfaceByDeCasteljau(
            in BezierSurface<Vector2> surface,
            Vector2 uv
        )
        {
            var degreeU = surface.DegreeU;
            var controlPoints = surface.ControlPoints;

            // Temp control points for the resulting curve in U direction
            var temp = new ControlPoint<Vector2>[degreeU + 1];

            for (int i = 0; i <= degreeU; i++)
            {
                // controlPoints[i] contains the row of points defining a curve in V direction.
                // We evaluate this curve at v (uv.Y).
                var curve = new BezierCurve<Vector2>(controlPoints[i]);
                var point = BezierCurve2DHelper.GetPointOnCurveByDeCasteljau(curve, uv.Y);
                temp[i] = new ControlPoint<Vector2>(point);
            }

            // evaluating the curve in U direction at u (uv.X)
            var bezierCurve = new BezierCurve<Vector2>(temp);
            return BezierCurve2DHelper.GetPointOnCurveByDeCasteljau(bezierCurve, uv.X);
        }
    }
}
