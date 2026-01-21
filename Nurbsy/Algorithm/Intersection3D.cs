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

namespace Nurbsy.Algorithm
{
    public static class Intersection3D
    {
        public static CurveCurveIntersectionType ComputeRays(
            Vector3 point0,
            Vector3 vector0,
            Vector3 point1,
            Vector3 vector1,
            out double param0,
            out double param1,
            out Vector3 intersectPoint
        )
        {
            if (MathUtils.IsAlmostEqualTo(vector0, vector1))
            {
                if (MathUtils.IsAlmostEqualTo(point0, point1))
                {
                    intersectPoint = point0;
                    param0 = param1 = 0;
                    return CurveCurveIntersectionType.Intersecting;
                }
            }

            Validate.Argument(
                MathUtils.IsZero(vector0),
                nameof(vector0),
                "Vector0 must not be zero vector."
            );
            Validate.Argument(
                MathUtils.IsZero(vector1),
                nameof(vector1),
                "Vector1 must not be zero vector."
            );

            Vector3 v0 = vector0;
            Vector3 v1 = vector1;

            Vector3 cross = Vector3.Cross(v0, v1);
            Vector3 diff = point1 - point0;
            Vector3 coinCross = Vector3.Cross(diff, v1);

            if (MathUtils.IsZero(cross))
            {
                if (MathUtils.IsZero(coinCross))
                {
                    intersectPoint = Vector3.Zero;
                    param0 = param1 = 0;
                    return CurveCurveIntersectionType.Coincident;
                }
                else
                {
                    intersectPoint = Vector3.Zero;
                    param0 = param1 = 0;
                    return CurveCurveIntersectionType.Parallel;
                }
            }

            double squareLength = cross.LengthSquared();

            Vector3 pd1Cross = Vector3.Cross(diff, v1);
            double pd1Dot = Vector3.Dot(pd1Cross, cross);
            param0 = pd1Dot / squareLength;

            Vector3 pd2Cross = Vector3.Cross(diff, v0);
            double pd2Dot = Vector3.Dot(pd2Cross, cross);
            param1 = pd2Dot / squareLength;

            Vector3 rayP0 = point0 + vector0 * (float)param0;
            Vector3 rayP1 = point1 + vector1 * (float)param1;

            if (MathUtils.IsAlmostEqualTo(rayP0, rayP1))
            {
                intersectPoint = rayP0;
                return CurveCurveIntersectionType.Intersecting;
            }

            intersectPoint = Vector3.Zero;
            return CurveCurveIntersectionType.Skew;
        }

        public static LinePlaneIntersectionType ComputeLineAndPlane(
            Vector3 normal,
            Vector3 pointOnPlane,
            Vector3 pointOnLine,
            Vector3 lineDirection,
            out Vector3 intersectPoint
        )
        {
            Vector3 planeNormal = normal;
            planeNormal.Normalize();
            Vector3 lineDirectionNormal = lineDirection;
            lineDirectionNormal.Normalize();
            Vector3 P2L = pointOnLine - pointOnPlane;
            Vector3 P2LNormal = P2L;
            if (!MathUtils.IsZero(P2LNormal))
                P2LNormal.Normalize();

            double dot = Vector3.Dot(P2LNormal, planeNormal);
            if (MathUtils.IsAlmostEqualTo(dot, 0.0))
            {
                intersectPoint = pointOnLine;
                return LinePlaneIntersectionType.On;
            }

            // Angle between plane normal and line direction
            double dotNormalDir = Vector3.Dot(planeNormal, lineDirectionNormal);
            // Clamp for acos
            dotNormalDir = Math.Clamp(dotNormalDir, -1.0, 1.0);
            double angle = Math.Acos(dotNormalDir);

            if (MathUtils.IsAlmostEqualTo(angle, Constants.Pi / 2))
            {
                intersectPoint = Vector3.Zero;
                return LinePlaneIntersectionType.Parallel;
            }

            double d = -Vector3.Dot(P2L, planeNormal) / Vector3.Dot(lineDirection, planeNormal);
            intersectPoint = (float)d * lineDirectionNormal + pointOnLine;
            return LinePlaneIntersectionType.Intersecting;
        }
    }
}
