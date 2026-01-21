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
    public static class Interpolation
    {
        public static double GetTotalChordLength<T>(IReadOnlyList<T> throughPoints)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                return Interpolation2D.GetTotalChordLength(input);
            }
            if (typeof(T) == typeof(Vector3))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                return Interpolation3D.GetTotalChordLength(input);
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static double[] GetChordParameterization<T>(IReadOnlyList<T> throughPoints)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                return Interpolation2D.GetChordParameterization(input);
            }
            if (typeof(T) == typeof(Vector3))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                return Interpolation3D.GetChordParameterization(input);
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static double GetCentripetalLength<T>(IReadOnlyList<T> throughPoints)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                return Interpolation2D.GetCentripetalLength(input);
            }
            if (typeof(T) == typeof(Vector3))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                return Interpolation3D.GetCentripetalLength(input);
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static double[] GetCentripetalParameterization<T>(IReadOnlyList<T> throughPoints)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                return Interpolation2D.GetCentripetalParameterization(input);
            }
            if (typeof(T) == typeof(Vector3))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                return Interpolation3D.GetCentripetalParameterization(input);
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        // This method relies on double[] and is generic in logic, keeping it here as shared
        public static double[] AverageKnotVector(int degree, IReadOnlyList<double> paramsList)
        {
            int size = paramsList.Count;
            int n = size - 1;
            int m = n + degree + 1;

            var knotVector = new double[m + 1];
            // C++: for (int i = m - degree; i <= m; i++) knotVector[i] = 1.0;
            // index m is included in size m+1
            // Careful about array size, Nurbs usually m = n + p + 1. size is m+1?
            // "vector<double> knotVector(m + 1, 0.0);"
            // It seems C++ implementation uses m+1 size. Let's follow.

            for (int i = m - degree; i <= m; i++)
            {
                knotVector[i] = 1.0;
            }

            for (int j = 1; j <= n - degree; j++)
            {
                double sum = 0.0;
                for (int i = j; i <= j + degree - 1; i++)
                {
                    sum += paramsList[i];
                }
                knotVector[j + degree] = (1.0 / degree) * sum;
            }
            return knotVector;
        }

        public static double[] ComputeKnotVector(
            int degree,
            int controlPointsCount,
            IReadOnlyList<double> paramsList
        )
        {
            int m = paramsList.Count - 1;
            int n = controlPointsCount - 1;
            int nn = n + degree + 2; // C++: nn = n + degree + 2, vector(nn)

            var knotVector = new double[nn];
            double d = (double)(m + 1) / (double)(n - degree + 1);

            for (int j = 1; j <= n - degree; j++)
            {
                int i = (int)Math.Floor(j * d);
                double alpha = (j * d) - i;
                // Check bounds for paramsList[i] and paramsList[i-1]
                // i is 1-based index from logical algo, but array 0-based?
                // C++: params[i - 1], params[i].
                knotVector[degree + j] =
                    (1.0 - alpha) * paramsList[i - 1] + (alpha * paramsList[i]);
            }

            for (int i = 0; i < nn; i++)
            {
                if (i <= degree)
                {
                    knotVector[i] = 0.0;
                }
                else if (i >= nn - 1 - degree)
                {
                    knotVector[i] = 1.0;
                }
            }
            return knotVector;
        }

        public static bool TryComputeTangents<T>(IReadOnlyList<T> throughPoints, out T[] tangents)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                if (Interpolation2D.TryComputeTangents(input, out var res))
                {
                    tangents = Unsafe.As<Vector2[], T[]>(ref res);
                    return true;
                }
                tangents = null;
                return false;
            }
            if (typeof(T) == typeof(Vector3))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                if (Interpolation3D.TryComputeTangents(input, out var res))
                {
                    tangents = Unsafe.As<Vector3[], T[]>(ref res);
                    return true;
                }
                tangents = null;
                return false;
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static T[] ComputeTangents<T>(IReadOnlyList<T> throughPoints)
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector2>>(ref throughPoints);
                var res = Interpolation2D.ComputeTangents(input);
                return Unsafe.As<Vector2[], T[]>(ref res);
            }
            if (typeof(T) == typeof(Vector3))
            {
                var input = Unsafe.As<IReadOnlyList<T>, IReadOnlyList<Vector3>>(ref throughPoints);
                var res = Interpolation3D.ComputeTangents(input);
                return Unsafe.As<Vector3[], T[]>(ref res);
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static bool ComputerWeightForRationalQuadraticInterpolation<T>(
            T startPoint,
            T middleControlPoint,
            T endPoint,
            out double weight
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                return Interpolation2D.ComputerWeightForRationalQuadraticInterpolation(
                    Unsafe.As<T, Vector2>(ref startPoint),
                    Unsafe.As<T, Vector2>(ref middleControlPoint),
                    Unsafe.As<T, Vector2>(ref endPoint),
                    out weight
                );
            }
            if (typeof(T) == typeof(Vector3))
            {
                return Interpolation3D.ComputerWeightForRationalQuadraticInterpolation(
                    Unsafe.As<T, Vector3>(ref startPoint),
                    Unsafe.As<T, Vector3>(ref middleControlPoint),
                    Unsafe.As<T, Vector3>(ref endPoint),
                    out weight
                );
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public static bool GetSurfaceMeshParameterization<T>(
            IReadOnlyList<IReadOnlyList<T>> throughPoints,
            out double[] paramsU,
            out double[] paramsV
        )
            where T : struct
        {
            if (typeof(T) == typeof(Vector2))
            {
                var input = Unsafe.As<
                    IReadOnlyList<IReadOnlyList<T>>,
                    IReadOnlyList<IReadOnlyList<Vector2>>
                >(ref throughPoints);
                return Interpolation2D.GetSurfaceMeshParameterization(
                    input,
                    out paramsU,
                    out paramsV
                );
            }
            if (typeof(T) == typeof(Vector3))
            {
                var input = Unsafe.As<
                    IReadOnlyList<IReadOnlyList<T>>,
                    IReadOnlyList<IReadOnlyList<Vector3>>
                >(ref throughPoints);
                return Interpolation3D.GetSurfaceMeshParameterization(
                    input,
                    out paramsU,
                    out paramsV
                );
            }
            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
