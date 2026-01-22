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
using Nurbsy.Helpers;

namespace Nurbsy
{
    public readonly record struct BezierCurve<T>
        where T : struct
    {
        /// <summary>
        /// The degree of the Bezier curve, calculated as (ControlPoints.Count - 1).
        /// </summary>
        public int Degree { get; }

        /// <summary>
        /// The control points defining the curve.
        /// </summary>
        public IReadOnlyList<ControlPoint<T>> ControlPoints { get; }

        /// <summary>
        /// Creates a Bezier curve from weighted control points.
        /// </summary>
        public BezierCurve(int degree, IReadOnlyList<ControlPoint<T>> controlPoints)
        {
            Degree = degree;
            ControlPoints = controlPoints;
            Check();
        }

        /// <summary>
        /// Creates a Bezier curve from unweighted points (weights set to 1.0).
        /// </summary>
        public BezierCurve(int degree, IReadOnlyList<T> controlPoints)
        {
            Degree = degree;
            ControlPoints = controlPoints.Select(cp => new ControlPoint<T>(cp)).ToImmutableArray();
            Check();
        }

        public void Check()
        {
            Validate.Argument(Degree > 0, nameof(Degree), "Degree must be greater than zero.");
            Validate.Argument(
                ControlPoints != null && ControlPoints.Count > 0,
                nameof(ControlPoints),
                "ControlPoints must contain at least one point."
            );
            Validate.Argument(
                Validate.IsValidBezier(Degree, ControlPoints.Count),
                nameof(ControlPoints),
                "ControlPoints count equals degree plus one."
            );
        }

        public static T GetPointOnQuadraticArc(
            ControlPoint<T> startPoint,
            ControlPoint<T> middlePoint,
            ControlPoint<T> endPoint,
            double paramT
        )
        {
            return BezierCurveHelper.GetPointOnQuadraticArc(
                startPoint,
                middlePoint,
                endPoint,
                paramT
            );
        }

        public static IReadOnlyList<ControlPoint<T>>? ComputeMiddleControlPointsOnQuadraticCurve(
            T startPoint,
            T startTangent,
            T endPoint,
            T endTangent
        )
        {
            return BezierCurveHelper.ComputeMiddleControlPointsOnQuadraticCurve(
                startPoint,
                startTangent,
                endPoint,
                endTangent
            );
        }

        public T GetPointOnCurveByBernstein(double paramT)
        {
            Validate.Range(paramT, 0.0, 1.0, nameof(paramT));
            return BezierCurveHelper.GetPointOnCurveByBernstein(this, paramT);
        }

        public T GetPointOnCurveByDeCasteljau(double paramT)
        {
            Validate.Range(paramT, 0.0, 1.0, nameof(paramT));
            return BezierCurveHelper.GetPointOnCurveByDeCasteljau(this, paramT);
        }
    }
}
