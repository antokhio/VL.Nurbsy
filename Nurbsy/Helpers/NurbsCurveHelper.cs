using System.Runtime.CompilerServices;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    public static class NurbsCurveHelper
    {
        public static NurbsCurve<T> CreateLine<T>(T start, T end)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var s = Unsafe.As<T, Vector2>(ref start);
                var e = Unsafe.As<T, Vector2>(ref end);
                var result = NurbsCurve2DHelper.CreateLine(s, e);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                var s = Unsafe.As<T, Vector3>(ref start);
                var e = Unsafe.As<T, Vector3>(ref end);
                var result = NurbsCurve3DHelper.CreateLine(s, e);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static NurbsCurve<T> CreateCubicHermite<T>(
            IReadOnlyList<T> throughPoints,
            IReadOnlyList<T> tangents
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var inputPoints = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(
                    ref throughPoints
                );
                var inputTangents = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(
                    ref tangents
                );
                var result = NurbsCurve2DHelper.CreateCubicHermite(inputPoints, inputTangents);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                var inputPoints = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(
                    ref throughPoints
                );
                var inputTangents = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(
                    ref tangents
                );
                var result = NurbsCurve3DHelper.CreateCubicHermite(inputPoints, inputTangents);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static NurbsCurve<T> CreateArc<T>(
            T center,
            T xAxis,
            T yAxis,
            double startRad,
            double endRad,
            double xRadius,
            double yRadius
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var c = Unsafe.As<T, Vector2>(ref center);
                var x = Unsafe.As<T, Vector2>(ref xAxis);
                var y = Unsafe.As<T, Vector2>(ref yAxis);
                var result = NurbsCurve2DHelper.CreateArc(
                    c,
                    x,
                    y,
                    startRad,
                    endRad,
                    xRadius,
                    yRadius
                );
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                var c = Unsafe.As<T, Vector3>(ref center);
                var x = Unsafe.As<T, Vector3>(ref xAxis);
                var y = Unsafe.As<T, Vector3>(ref yAxis);
                var result = NurbsCurve3DHelper.CreateArc(
                    c,
                    x,
                    y,
                    startRad,
                    endRad,
                    xRadius,
                    yRadius
                );
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static bool CreateOneConicArc<T>(
            T start,
            T startTangent,
            T end,
            T endTangent,
            T pointOnConic,
            out T projectPoint,
            out double projectPointWeight
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var s = Unsafe.As<T, Vector2>(ref start);
                var st = Unsafe.As<T, Vector2>(ref startTangent);
                var e = Unsafe.As<T, Vector2>(ref end);
                var et = Unsafe.As<T, Vector2>(ref endTangent);
                var pc = Unsafe.As<T, Vector2>(ref pointOnConic);

                bool result = NurbsCurve2DHelper.CreateOneConicArc(
                    s,
                    st,
                    e,
                    et,
                    pc,
                    out var pp,
                    out projectPointWeight
                );
                projectPoint = Unsafe.As<Vector2, T>(ref pp);
                return result;
            }

            if (typeof(T) == typeof(Vector3))
            {
                var s = Unsafe.As<T, Vector3>(ref start);
                var st = Unsafe.As<T, Vector3>(ref startTangent);
                var e = Unsafe.As<T, Vector3>(ref end);
                var et = Unsafe.As<T, Vector3>(ref endTangent);
                var pc = Unsafe.As<T, Vector3>(ref pointOnConic);

                bool result = NurbsCurve3DHelper.CreateOneConicArc(
                    s,
                    st,
                    e,
                    et,
                    pc,
                    out var pp,
                    out projectPointWeight
                );
                projectPoint = Unsafe.As<Vector3, T>(ref pp);
                return result;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static bool CreateOpenConic<T>(
            T start,
            T startTangent,
            T end,
            T endTangent,
            T pointOnConic,
            out NurbsCurve<T> curve
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var s = Unsafe.As<T, Vector2>(ref start);
                var st = Unsafe.As<T, Vector2>(ref startTangent);
                var e = Unsafe.As<T, Vector2>(ref end);
                var et = Unsafe.As<T, Vector2>(ref endTangent);
                var pc = Unsafe.As<T, Vector2>(ref pointOnConic);

                if (NurbsCurve2DHelper.CreateOpenConic(s, st, e, et, pc, out var res2))
                {
                    curve = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                curve = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                var s = Unsafe.As<T, Vector3>(ref start);
                var st = Unsafe.As<T, Vector3>(ref startTangent);
                var e = Unsafe.As<T, Vector3>(ref end);
                var et = Unsafe.As<T, Vector3>(ref endTangent);
                var pc = Unsafe.As<T, Vector3>(ref pointOnConic);

                if (NurbsCurve3DHelper.CreateOpenConic(s, st, e, et, pc, out var res3))
                {
                    curve = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                curve = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static NurbsCurve<T> GlobalInterpolation<T>(
            int degree,
            IReadOnlyList<T> throughPoints,
            IReadOnlyList<double> parameters = null
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var points = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                var result = NurbsCurve2DHelper.GlobalInterpolation(degree, points, parameters);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                var points = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                var result = NurbsCurve3DHelper.GlobalInterpolation(degree, points, parameters);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static NurbsCurve<T> GlobalInterpolation<T>(
            int degree,
            IReadOnlyList<T> throughPoints,
            IReadOnlyList<T> tangents,
            double tangentFactor = 1.0
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var points = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                var tang = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref tangents);
                var result = NurbsCurve2DHelper.GlobalInterpolation(
                    degree,
                    points,
                    tang,
                    tangentFactor
                );
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                var points = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                var tang = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref tangents);
                var result = NurbsCurve3DHelper.GlobalInterpolation(
                    degree,
                    points,
                    tang,
                    tangentFactor
                );
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static bool WeightedAndContrainedLeastSquaresApproximation<T>(
            int degree,
            IReadOnlyList<T> throughPoints,
            IReadOnlyList<double> throughPointWeights,
            IReadOnlyList<T> tangents,
            IReadOnlyList<int> tangentIndices,
            IReadOnlyList<double> tangentWeights,
            int controlPointsCount,
            out NurbsCurve<T> curve
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var tp = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                var tg = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref tangents);

                if (
                    NurbsCurve2DHelper.WeightedAndContrainedLeastSquaresApproximation(
                        degree,
                        tp,
                        throughPointWeights,
                        tg,
                        tangentIndices,
                        tangentWeights,
                        controlPointsCount,
                        out var res2
                    )
                )
                {
                    curve = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                curve = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                var tp = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                var tg = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref tangents);

                if (
                    NurbsCurve3DHelper.WeightedAndContrainedLeastSquaresApproximation(
                        degree,
                        tp,
                        throughPointWeights,
                        tg,
                        tangentIndices,
                        tangentWeights,
                        controlPointsCount,
                        out var res3
                    )
                )
                {
                    curve = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                curve = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static bool GlobalApproximationByErrorBound<T>(
            int degree,
            IReadOnlyList<T> throughPoints,
            double maxError,
            out NurbsCurve<T> curve
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var points = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                if (
                    NurbsCurve2DHelper.GlobalApproximationByErrorBound(
                        degree,
                        points,
                        maxError,
                        out var res2
                    )
                )
                {
                    curve = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                curve = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                var points = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                if (
                    NurbsCurve3DHelper.GlobalApproximationByErrorBound(
                        degree,
                        points,
                        maxError,
                        out var res3
                    )
                )
                {
                    curve = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                curve = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static void SplitArc<T>(
            T start,
            T projectPoint,
            double projectPointWeight,
            T end,
            out T insertPointAtStartSide,
            out T splitPoint,
            out T insertPointAtEndSide,
            out double insertWeight
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var s = Unsafe.As<T, Vector2>(ref start);
                var pp = Unsafe.As<T, Vector2>(ref projectPoint);
                var e = Unsafe.As<T, Vector2>(ref end);

                NurbsCurve2DHelper.SplitArc(
                    s,
                    pp,
                    projectPointWeight,
                    e,
                    out var startSide,
                    out var sp,
                    out var endSide,
                    out insertWeight
                );

                insertPointAtStartSide = Unsafe.As<Vector2, T>(ref startSide);
                splitPoint = Unsafe.As<Vector2, T>(ref sp);
                insertPointAtEndSide = Unsafe.As<Vector2, T>(ref endSide);
                return;
            }

            if (typeof(T) == typeof(Vector3))
            {
                var s = Unsafe.As<T, Vector3>(ref start);
                var pp = Unsafe.As<T, Vector3>(ref projectPoint);
                var e = Unsafe.As<T, Vector3>(ref end);

                NurbsCurve3DHelper.SplitArc(
                    s,
                    pp,
                    projectPointWeight,
                    e,
                    out var startSide,
                    out var sp,
                    out var endSide,
                    out insertWeight
                );

                insertPointAtStartSide = Unsafe.As<Vector3, T>(ref startSide);
                splitPoint = Unsafe.As<Vector3, T>(ref sp);
                insertPointAtEndSide = Unsafe.As<Vector3, T>(ref endSide);
                return;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
