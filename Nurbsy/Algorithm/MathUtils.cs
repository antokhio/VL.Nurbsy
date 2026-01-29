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

using System.Diagnostics;
using System.Runtime.CompilerServices;
using Stride.Core.Mathematics;

namespace Nurbsy.Algorithm
{
    internal static class MathUtils
    {
        public const double Epsilon = 1e-9;

        public static bool IsAlmostEqualTo(double a, double b, double epsilon = Epsilon)
        {
            return Math.Abs(a - b) <= epsilon;
        }

        public static bool IsAlmostEqualTo(float a, float b, double epsilon = Epsilon)
        {
            return Math.Abs(a - b) <= epsilon;
        }

        public static bool IsAlmostEqualTo(Vector2 a, Vector2 b, double epsilon = Epsilon)
        {
            return Vector2.DistanceSquared(a, b) <= (float)(epsilon * epsilon);
        }

        public static bool IsAlmostEqualTo(Vector3 a, Vector3 b, double epsilon = Epsilon)
        {
            return Vector3.DistanceSquared(a, b) <= (float)(epsilon * epsilon);
        }

        public static bool IsGreaterThanOrEqual(double a, double b, double epsilon = Epsilon)
        {
            return a >= b - epsilon;
        }

        public static bool IsLessThanOrEqual(double a, double b, double epsilon = Epsilon)
        {
            return a <= b + epsilon;
        }

        public static bool IsLessThan(double a, double b, double epsilon = Epsilon)
        {
            return a < b - epsilon;
        }

        public static bool IsGreaterThan(double a, double b, double epsilon = Epsilon)
        {
            return a > b + epsilon;
        }

        public static bool IsZero(double a, double epsilon = Epsilon)
        {
            return Math.Abs(a) <= epsilon;
        }

        public static bool IsZero(Vector2 v, double epsilon = Epsilon)
        {
            return v.LengthSquared() <= (float)(epsilon * epsilon);
        }

        public static bool IsZero(Vector3 v, double epsilon = Epsilon)
        {
            return v.LengthSquared() <= (float)(epsilon * epsilon);
        }

        // Gaussian elimination solver for Ax = B
        public static double[][] SolveLinearSystem(double[][] A, double[][] B)
        {
            int n = A.Length;
            int m = A[0].Length;
            int p = B[0].Length;

            Debug.Assert(n == m, "Matrix A must be square.");
            Debug.Assert(A.Length == B.Length, "Matrix A and B must have same number of rows.");

            // Deep copy A and B to avoid modifying originals (if needed, or just work on clones)
            // Working with arrays of arrays
            var ACopy = new double[n][];
            for (int i = 0; i < n; i++)
            {
                ACopy[i] = new double[n];
                Array.Copy(A[i], ACopy[i], n);
            }

            var Result = new double[n][];
            for (int i = 0; i < n; i++)
            {
                Result[i] = new double[p];
                Array.Copy(B[i], Result[i], p);
            }

            // Forward elimination
            for (int i = 0; i < n; i++)
            {
                // Pivot
                int pivot = i;
                for (int j = i + 1; j < n; j++)
                {
                    if (Math.Abs(ACopy[j][i]) > Math.Abs(ACopy[pivot][i]))
                    {
                        pivot = j;
                    }
                }

                // Swap rows
                var tempA = ACopy[i];
                ACopy[i] = ACopy[pivot];
                ACopy[pivot] = tempA;

                var tempB = Result[i];
                Result[i] = Result[pivot];
                Result[pivot] = tempB;

                if (Math.Abs(ACopy[i][i]) < Epsilon)
                    throw new InvalidOperationException("Matrix is singular.");

                for (int j = i + 1; j < n; j++)
                {
                    double factor = ACopy[j][i] / ACopy[i][i];
                    for (int k = i; k < n; k++)
                    {
                        ACopy[j][k] -= factor * ACopy[i][k];
                    }
                    for (int k = 0; k < p; k++)
                    {
                        Result[j][k] -= factor * Result[i][k];
                    }
                }
            }

            // Backward substitution
            for (int i = n - 1; i >= 0; i--)
            {
                for (int k = 0; k < p; k++)
                {
                    double sum = 0;
                    for (int j = i + 1; j < n; j++)
                    {
                        sum += ACopy[i][j] * Result[j][k];
                    }
                    Result[i][k] = (Result[i][k] - sum) / ACopy[i][i];
                }
            }

            return Result;
        }

        public static double Binomial(int n, int k)
        {
            if (k < 0 || k > n)
                return 0;
            if (k == 0 || k == n)
                return 1;
            if (k > n / 2)
                k = n - k;

            double res = 1;
            for (int i = 1; i <= k; ++i)
            {
                res = res * (n - i + 1) / i;
            }
            return res;
        }

        public static double[][] Transpose(double[][] A)
        {
            int rows = A.Length;
            int cols = A[0].Length;
            var result = new double[cols][];
            for (int i = 0; i < cols; i++)
            {
                result[i] = new double[rows];
                for (int j = 0; j < rows; j++)
                {
                    result[i][j] = A[j][i];
                }
            }
            return result;
        }

        public static double[][] MatrixMultiply(double[][] A, double[][] B)
        {
            int rowsA = A.Length;
            int colsA = A[0].Length;
            int rowsB = B.Length;
            int colsB = B[0].Length;

            if (colsA != rowsB)
                throw new ArgumentException("Matrix dimensions do not match for multiplication.");

            var result = new double[rowsA][];
            for (int i = 0; i < rowsA; i++)
            {
                result[i] = new double[colsB];
                for (int j = 0; j < colsB; j++)
                {
                    double sum = 0.0;
                    for (int k = 0; k < colsA; k++)
                    {
                        sum += A[i][k] * B[k][j];
                    }
                    result[i][j] = sum;
                }
            }
            return result;
        }

        public static double[][] Identity(int n)
        {
            var I = new double[n][];
            for (int i = 0; i < n; i++)
            {
                I[i] = new double[n];
                I[i][i] = 1.0;
            }
            return I;
        }

        public static bool MakeInverse(double[][] A, out double[][] inverse)
        {
            inverse = null;
            try
            {
                int n = A.Length;
                var identity = Identity(n);
                // Solving Ax = I gives x = A^-1
                inverse = SolveLinearSystem(A, identity);
                return true;
            }
            catch
            {
                return false;
            }
        }

        public static double[][] MakeDiagonal(IReadOnlyList<double> diag)
        {
            int n = diag.Count;
            var M = new double[n][];
            for (int i = 0; i < n; i++)
            {
                M[i] = new double[n];
                M[i][i] = diag[i];
            }
            return M;
        }

        public static IReadOnlyList<IReadOnlyList<T>> Transpose<T>(
            IReadOnlyList<IReadOnlyList<T>> input
        )
        {
            int rows = input.Count;
            int cols = input[0].Count;
            var result = new T[cols][];
            for (int i = 0; i < cols; i++)
            {
                result[i] = new T[rows];
                for (int j = 0; j < rows; j++)
                {
                    result[i][j] = input[j][i];
                }
            }
            return result;
        }

        internal static double Distance<T>(T a, T b)
        {
            if (typeof(T) == typeof(Vector3))
            {
                var a3 = Unsafe.As<T, Vector3>(ref a);
                var b3 = Unsafe.As<T, Vector3>(ref b);

                return Vector3.Distance(a3, b3);
            }
            if (typeof(T) == typeof(Vector2))
            {
                var a2 = Unsafe.As<T, Vector2>(ref a);
                var b2 = Unsafe.As<T, Vector2>(ref b);

                return Vector2.Distance(a2, b2);
            }
            return 0.0;
        }
    }
}
