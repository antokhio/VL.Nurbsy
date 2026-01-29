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

using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    public static class NurbsSurfaceHelper
    {
        /// <summary>
        /// The NURBS Book 2nd Edition Page337
        /// Create a ruled surface.
        /// </summary>
        public static NurbsSurface<Vector3> CreateRulledSurface(
            NurbsCurve<Vector3> curve1,
            NurbsCurve<Vector3> curve2
        )
        {
            var k1 = curve1.Knots;
            var k2 = curve2.Knots;

            // 1. Validate parameter ranges
            bool rangesMatch =
                MathUtils.IsAlmostEqualTo(k1[0], k2[0])
                && MathUtils.IsAlmostEqualTo(k1[k1.Count - 1], k2[k2.Count - 1]);

            Validate.Argument(
                rangesMatch,
                "curve1 & curve2",
                "Ensure that the two curves are defined on the same parameter range."
            );

            // 2. Elevate degrees to match the maximum degree
            int degree1 = curve1.Degree;
            int degree2 = curve2.Degree;
            int maxDegree = Math.Max(degree1, degree2);

            if (degree1 < maxDegree)
                curve1 = curve1.ElevateDegree(maxDegree - degree1);

            if (degree2 < maxDegree)
                curve2 = curve2.ElevateDegree(maxDegree - degree2);

            // 3. Unify knot vectors (merged sets)
            KnotsUtils.GetInsertionKnots(
                curve1.Knots,
                curve2.Knots,
                out var insert1,
                out var insert2
            );

            if (insert1.Count > 0)
                curve1 = curve1.RefineKnotVector(insert1);

            if (insert2.Count > 0)
                curve2 = curve2.RefineKnotVector(insert2);

            // 4. Construct the ruled surface
            // The U direction follows the curves' degree and knots.
            // The V direction is linear (Degree 1) connecting the two curves.
            int finalDegreeU = maxDegree;
            int finalDegreeV = 1;

            var knotsU = curve1.Knots;
            // Standard clamped knot vector for degree 1: {0, 0, 1, 1}
            // (Assuming V parameterization 0 to 1)
            var knotsV = new double[] { 0.0, 0.0, 1.0, 1.0 };

            int uCount = curve1.ControlPoints.Count;
            var finalControlPoints = new ControlPoint<Vector3>[uCount][];

            for (int i = 0; i < uCount; i++)
            {
                // Surface[u][v] -> Row i corresponds to u_i, Columns 0,1 correspond to curve1, curve2
                finalControlPoints[i] = new ControlPoint<Vector3>[]
                {
                    curve1.ControlPoints[i],
                    curve2.ControlPoints[i],
                };
            }

            return new NurbsSurface<Vector3>(
                finalDegreeU,
                finalDegreeV,
                finalControlPoints,
                knotsU,
                knotsV
            );
        }

        /// <summary>
        /// Create a loft surface from a list of sections (NURBS curves).
        /// </summary>
        public static NurbsSurface<T> CreateLoftSurface<T>(
            IReadOnlyList<NurbsCurve<T>> sections,
            int? degreeV = null,
            IReadOnlyList<double> customTrajectoryParameters = null
        )
            where T : struct
        {
            Validate.Argument(
                sections != null && sections.Count > 0,
                nameof(sections),
                "Sections must not be empty."
            );

            // 1. Unify curves (Degree and Knots)
            int maxDegree = 0;
            foreach (var curve in sections)
            {
                if (curve.Degree > maxDegree)
                    maxDegree = curve.Degree;
            }

            var unifiedCurves = new List<NurbsCurve<T>>();
            foreach (var curve in sections)
            {
                var current = curve;
                if (current.Degree < maxDegree)
                {
                    current = current.ElevateDegree(maxDegree - current.Degree);
                }
                unifiedCurves.Add(current);
            }

            var commonKnots = new List<double>(unifiedCurves[0].Knots);
            for (int i = 1; i < unifiedCurves.Count; i++)
            {
                KnotsUtils.GetInsertionKnots(
                    commonKnots,
                    unifiedCurves[i].Knots,
                    out var insert1,
                    out _
                );
                if (insert1.Count > 0)
                {
                    commonKnots.AddRange(insert1);
                    commonKnots.Sort();
                }
            }

            for (int i = 0; i < unifiedCurves.Count; i++)
            {
                KnotsUtils.GetInsertionKnots(
                    unifiedCurves[i].Knots,
                    commonKnots,
                    out var insert,
                    out _
                );
                if (insert.Count > 0)
                {
                    unifiedCurves[i] = unifiedCurves[i].RefineKnotVector(insert);
                }
            }

            // 2. Compute parameters vl
            int finalDegreeV = degreeV ?? maxDegree;
            IReadOnlyList<double> vl;

            if (customTrajectoryParameters != null)
            {
                vl = customTrajectoryParameters;
            }
            else
            {
                int n = unifiedCurves[0].ControlPoints.Count - 1; // CP count - 1
                int kSize = unifiedCurves.Count;
                var computedVl = new double[kSize];
                computedVl[0] = 0.0;
                computedVl[kSize - 1] = 1.0;

                // Precompute lengths of each column track
                var columnLengths = new double[n + 1];
                for (int i = 0; i <= n; i++)
                {
                    var colPoints = new T[kSize];
                    for (int k = 0; k < kSize; k++)
                    {
                        colPoints[k] = unifiedCurves[k].ControlPoints[i].Value;
                    }
                    columnLengths[i] = Interpolation.GetTotalChordLength(
                        (IReadOnlyList<T>)colPoints
                    );
                }

                for (int k = 1; k < kSize - 1; k++)
                {
                    double sum = 0.0;
                    for (int i = 0; i <= n; i++)
                    {
                        double delta = MathUtils.Distance(
                            unifiedCurves[k].ControlPoints[i].Value,
                            unifiedCurves[k - 1].ControlPoints[i].Value
                        );
                        if (columnLengths[i] > MathUtils.Epsilon)
                        {
                            sum += delta / columnLengths[i];
                        }
                    }
                    computedVl[k] = computedVl[k - 1] + (1.0 / (n + 1)) * sum;
                }
                vl = computedVl;
            }

            var knotsV = Interpolation.AverageKnotVector(finalDegreeV, vl);

            // 3. Interpolate columns
            int rows = unifiedCurves[0].ControlPoints.Count;
            int cols = unifiedCurves.Count;

            var resultControlPoints = new ControlPoint<T>[rows][];

            for (int i = 0; i < rows; i++)
            {
                var columnPoints = new T[cols];
                for (int k = 0; k < cols; k++)
                {
                    columnPoints[k] = unifiedCurves[k].ControlPoints[i].Value;
                }

                var curve = NurbsCurveHelper.GlobalInterpolation(finalDegreeV, columnPoints, vl);

                var resultRow = new ControlPoint<T>[cols];
                for (int k = 0; k < cols; k++)
                {
                    T val = curve.ControlPoints[k].Value;
                    double weight = unifiedCurves[k].ControlPoints[i].Weight;
                    resultRow[k] = new ControlPoint<T>(val, weight);
                }
                resultControlPoints[i] = resultRow;
            }

            return new NurbsSurface<T>(
                maxDegree,
                finalDegreeV,
                resultControlPoints,
                commonKnots,
                knotsV
            );
        }
    }
}
