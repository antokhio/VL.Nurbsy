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
    public enum CurveCurveIntersectionType
    {
        Intersecting = 0,
        Parallel = 1,
        Coincident = 2,
        Skew = 3,
    }

    public enum LinePlaneIntersectionType
    {
        Intersecting = 0,
        Parallel = 1,
        On = 2,
    }

    public enum CurveNormal
    {
        Normal = 0,
        Binormal = 1,
    }

    public enum SurfaceDirection
    {
        All = 0,
        UDirection = 1,
        VDirection = 2,
    }

    /// <summary>
    /// Types of surface curvature.
    /// </summary>
    public enum SurfaceCurvature
    {
        /// <summary>
        /// Gaussian curvature (K = k1 * k2)
        /// </summary>
        Gauss,

        /// <summary>
        /// Mean curvature (H = (k1 + k2) / 2)
        /// </summary>
        Mean,

        /// <summary>
        /// Maximum principal curvature (k1)
        /// </summary>
        Maximum,

        /// <summary>
        /// Minimum principal curvature (k2)
        /// </summary>
        Minimum,

        /// <summary>
        /// Absolute curvature (|k1| + |k2|)
        /// </summary>
        Abs,

        /// <summary>
        /// Root mean square curvature (sqrt(k1² + k2²))
        /// </summary>
        Rms,
    }

    public enum IntegratorType
    {
        Simpson = 0,
        GaussLegendre = 1,
        Chebyshev = 2,
    }

    public enum OffsetType
    {
        TillerAndHanson = 0,
        PieglAndTiller = 1,
    }

    /// <summary>
    /// Specifies the algorithm used to evaluate a point on a Bézier curve at parameter t.
    /// </summary>
    public enum BezierEvaluation
    {
        /// <summary>
        /// Evaluates the curve using the explicit Bernstein polynomial definition (Weighted Sum).
        /// Faster execution, may lose precision with very high-degree curves due to large binomial coefficients.
        /// </summary>
        Bernstein,

        /// <summary>
        /// Evaluates the curve using De Casteljau's algorithm (Recursive Linear Interpolation).
        /// Slower execution. Numerically most stable method.
        /// </summary>
        DeCasteljau,
    }
}
