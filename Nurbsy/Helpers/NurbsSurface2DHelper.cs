using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    public static class NurbsSurface2DHelper
    {
        /// <inheritdoc cref="NurbsSurface{T}.GetPointOnSurface(Vector2)"/>
        public static Vector2 GetPointOnSurface(in NurbsSurface<Vector2> surface, Vector2 uv)
        {
            var degreeU = surface.DegreeU;
            var degreeV = surface.DegreeV;
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;
            var controlPoints = surface.ControlPoints;

            Validate.Range(uv.X, knotsU[0], knotsU[knotsU.Count - 1], nameof(uv.X));
            Validate.Range(uv.Y, knotsV[0], knotsV[knotsV.Count - 1], nameof(uv.Y));

            int uSpanIndex = Polynomials.GetKnotSpanIndex(degreeU, knotsU, uv.X);
            double[] Nu = Polynomials.BasisFunctions(uSpanIndex, degreeU, knotsU, uv.X);

            int vSpanIndex = Polynomials.GetKnotSpanIndex(degreeV, knotsV, uv.Y);
            double[] Nv = Polynomials.BasisFunctions(vSpanIndex, degreeV, knotsV, uv.Y);

            int uind = uSpanIndex - degreeU;

            // Use ControlPoint to accumulate (WeightedValue, Weight)
            var result = new ControlPoint<Vector2>(Vector2.Zero, 0.0);

            for (int l = 0; l <= degreeV; l++)
            {
                var temp = new ControlPoint<Vector2>(Vector2.Zero, 0.0);
                int vind = vSpanIndex - degreeV + l;

                for (int k = 0; k <= degreeU; k++)
                {
                    var cp = controlPoints[uind + k][vind];
                    double weight = cp.Weight;
                    double basisU = Nu[k];

                    temp.Value += cp.Value * (float)(weight * basisU);
                    temp.Weight += weight * basisU;
                }

                double basisV = Nv[l];
                result.Value += temp.Value * (float)basisV;
                result.Weight += temp.Weight * basisV;
            }

            return result.Weight != 0 ? result.Value / (float)result.Weight : Vector2.Zero;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetDerivatives(int, Vector2)"/>
        public static Vector2[][] ComputeDerivatives(
            in NurbsSurface<Vector2> surface,
            int derivative,
            Vector2 uv
        )
        {
            var degreeU = surface.DegreeU;
            var degreeV = surface.DegreeV;
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;
            var controlPoints = surface.ControlPoints;

            Validate.Argument(
                derivative > 0,
                nameof(derivative),
                "Derivative must be greater than zero."
            );
            Validate.Range(uv.X, knotsU[0], knotsU[knotsU.Count - 1], nameof(uv.X));
            Validate.Range(uv.Y, knotsV[0], knotsV[knotsV.Count - 1], nameof(uv.Y));

            // Initialize result array
            var derivatives = new Vector2[derivative + 1][];
            for (int i = 0; i <= derivative; i++)
            {
                derivatives[i] = new Vector2[derivative + 1];
            }

            int uSpanIndex = Polynomials.GetKnotSpanIndex(degreeU, knotsU, uv.X);
            double[][] Nu = Polynomials.BasisFunctionsDerivatives(
                uSpanIndex,
                degreeU,
                derivative,
                knotsU,
                uv.X
            );

            int vSpanIndex = Polynomials.GetKnotSpanIndex(degreeV, knotsV, uv.Y);
            double[][] Nv = Polynomials.BasisFunctionsDerivatives(
                vSpanIndex,
                degreeV,
                derivative,
                knotsV,
                uv.Y
            );

            int du = Math.Min(derivative, degreeU);
            int dv = Math.Min(derivative, degreeV);

            var temp = new Vector2[degreeV + 1];

            for (int k = 0; k <= du; k++)
            {
                for (int s = 0; s <= degreeV; s++)
                {
                    temp[s] = Vector2.Zero;
                    for (int r = 0; r <= degreeU; r++)
                    {
                        var cp = controlPoints[uSpanIndex - degreeU + r][vSpanIndex - degreeV + s];
                        temp[s] += cp.Value * (float)(cp.Weight * Nu[k][r]);
                    }
                }

                int dd = Math.Min(derivative - k, dv);
                for (int l = 0; l <= dd; l++)
                {
                    for (int s = 0; s <= degreeV; s++)
                    {
                        derivatives[k][l] += temp[s] * (float)Nv[l][s];
                    }
                }
            }

            return derivatives;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetFirstOrderDerivatives(Vector2)"/>
        public static Vector2[][] ComputeFirstOrderDerivatives(
            in NurbsSurface<Vector2> surface,
            Vector2 uv
        )
        {
            var degreeU = surface.DegreeU;
            var degreeV = surface.DegreeV;
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;
            var controlPoints = surface.ControlPoints;

            Validate.Range(uv.X, knotsU[0], knotsU[knotsU.Count - 1], nameof(uv.X));
            Validate.Range(uv.Y, knotsV[0], knotsV[knotsV.Count - 1], nameof(uv.Y));

            var derivatives = new Vector2[2][] { new Vector2[2], new Vector2[2] };

            int uSpanIndex = Polynomials.GetKnotSpanIndex(degreeU, knotsU, uv.X);
            double[][] Nu = Polynomials.BasisFunctionsFirstOrderDerivative(
                uSpanIndex,
                degreeU,
                knotsU,
                uv.X
            );

            int vSpanIndex = Polynomials.GetKnotSpanIndex(degreeV, knotsV, uv.Y);
            double[][] Nv = Polynomials.BasisFunctionsFirstOrderDerivative(
                vSpanIndex,
                degreeV,
                knotsV,
                uv.Y
            );

            int du = Math.Min(1, degreeU);
            int dv = Math.Min(1, degreeV);

            var temp = new Vector2[degreeV + 1];

            for (int k = 0; k <= du; k++)
            {
                for (int s = 0; s <= degreeV; s++)
                {
                    temp[s] = Vector2.Zero;
                    for (int r = 0; r <= degreeU; r++)
                    {
                        var cp = controlPoints[uSpanIndex - degreeU + r][vSpanIndex - degreeV + s];
                        temp[s] += cp.Value * (float)(cp.Weight * Nu[k][r]);
                    }
                }

                int dd = Math.Min(1, dv);
                for (int l = 0; l <= dd; l++)
                {
                    for (int s = 0; s <= degreeV; s++)
                    {
                        derivatives[k][l] += temp[s] * (float)Nv[l][s];
                    }
                }
            }

            return derivatives;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetControlPointsOfDerivatives(int, int, int, int, int)"/>
        public static Vector2[][][][] ComputeControlPointsOfDerivatives(
            in NurbsSurface<Vector2> surface,
            int derivative,
            int minSpanIndexU,
            int maxSpanIndexU,
            int minSpanIndexV,
            int maxSpanIndexV
        )
        {
            var degreeU = surface.DegreeU;
            var degreeV = surface.DegreeV;
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;
            var controlPoints = surface.ControlPoints;

            Validate.Argument(
                derivative > 0,
                nameof(derivative),
                "Derivative must be greater than zero."
            );
            Validate.Argument(
                minSpanIndexU >= 0 && minSpanIndexU <= maxSpanIndexU,
                nameof(minSpanIndexU),
                "Invalid span index range."
            );
            Validate.Argument(
                minSpanIndexV >= 0 && minSpanIndexV <= maxSpanIndexV,
                nameof(minSpanIndexV),
                "Invalid span index range."
            );

            int du = Math.Min(derivative, degreeU);
            int dv = Math.Min(derivative, degreeV);
            int rangeU = maxSpanIndexU - minSpanIndexU;
            int rangeV = maxSpanIndexV - minSpanIndexV;

            // Initialize PKL[k][l][i][j]
            var PKL = new Vector2[derivative + 1][][][];
            for (int k = 0; k <= derivative; k++)
            {
                PKL[k] = new Vector2[derivative + 1][][];
                for (int l = 0; l <= derivative; l++)
                {
                    PKL[k][l] = new Vector2[rangeU + 1][];
                    for (int i = 0; i <= rangeU; i++)
                    {
                        PKL[k][l][i] = new Vector2[rangeV + 1];
                    }
                }
            }

            // Compute derivatives in U direction for each V column
            for (int j = minSpanIndexV; j <= maxSpanIndexV; j++)
            {
                var points = new ControlPoint<Vector2>[controlPoints.Count];
                for (int i = 0; i < controlPoints.Count; i++)
                {
                    points[i] = controlPoints[i][j];
                }

                var tempCurve = new NurbsCurve<Vector2>(degreeU, points, knotsU);
                var temp = NurbsCurve2DHelper.ComputeControlPointsOfDerivatives(
                    in tempCurve,
                    du,
                    minSpanIndexU,
                    maxSpanIndexU
                );

                for (int k = 0; k <= du; k++)
                {
                    for (int i = 0; i <= rangeU - k; i++)
                    {
                        PKL[k][0][i][j - minSpanIndexV] = temp[k][i];
                    }
                }
            }

            // Compute derivatives in V direction
            for (int k = 0; k <= du; k++)
            {
                for (int i = 0; i <= rangeU - k; i++)
                {
                    int dd = Math.Min(derivative - k, dv);

                    var points = new ControlPoint<Vector2>[rangeV + 1];
                    for (int j = 0; j <= rangeV; j++)
                    {
                        points[j] = new ControlPoint<Vector2>(PKL[k][0][i][j], 1.0);
                    }

                    var tempKnots = new double[rangeV + degreeV + 2];
                    for (
                        int idx = 0;
                        idx < tempKnots.Length && minSpanIndexV + idx < knotsV.Count;
                        idx++
                    )
                    {
                        tempKnots[idx] = knotsV[minSpanIndexV + idx];
                    }

                    var tempCurve = new NurbsCurve<Vector2>(degreeV, points, tempKnots);
                    var temp = NurbsCurve2DHelper.ComputeControlPointsOfDerivatives(
                        in tempCurve,
                        dd,
                        0,
                        rangeV
                    );

                    for (int l = 1; l <= dd; l++)
                    {
                        for (int j = 0; j <= rangeV - l; j++)
                        {
                            PKL[k][l][i][j] = temp[l][j];
                        }
                    }
                }
            }

            return PKL;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetDerivativesByAllBasisFunctions(int, Vector2)"/>
        public static Vector2[][] ComputeDerivativesByAllBasisFunctions(
            in NurbsSurface<Vector2> surface,
            int derivative,
            Vector2 uv
        )
        {
            var degreeU = surface.DegreeU;
            var degreeV = surface.DegreeV;
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;

            Validate.Argument(
                derivative > 0,
                nameof(derivative),
                "Derivative must be greater than zero."
            );
            Validate.Range(uv.X, knotsU[0], knotsU[knotsU.Count - 1], nameof(uv.X));
            Validate.Range(uv.Y, knotsV[0], knotsV[knotsV.Count - 1], nameof(uv.Y));

            var SKL = new Vector2[derivative + 1][];
            for (int i = 0; i <= derivative; i++)
            {
                SKL[i] = new Vector2[derivative + 1];
            }

            int uSpanIndex = Polynomials.GetKnotSpanIndex(degreeU, knotsU, uv.X);
            int vSpanIndex = Polynomials.GetKnotSpanIndex(degreeV, knotsV, uv.Y);
            double[][] Nu = Polynomials.AllBasisFunctions(uSpanIndex, degreeU, knotsU, uv.X);
            double[][] Nv = Polynomials.AllBasisFunctions(vSpanIndex, degreeV, knotsV, uv.Y);

            var PKL = ComputeControlPointsOfDerivatives(
                in surface,
                derivative,
                uSpanIndex - degreeU,
                uSpanIndex,
                vSpanIndex - degreeV,
                vSpanIndex
            );

            int du = Math.Min(derivative, degreeU);
            int dv = Math.Min(derivative, degreeV);

            for (int k = 0; k <= du; k++)
            {
                int dd = Math.Min(derivative - k, dv);
                for (int l = 0; l <= dd; l++)
                {
                    SKL[k][l] = Vector2.Zero;
                    for (int i = 0; i <= degreeV - l; i++)
                    {
                        var temp = Vector2.Zero;
                        for (int j = 0; j <= degreeU - k; j++)
                        {
                            temp += PKL[k][l][j][i] * (float)Nu[j][degreeU - k];
                        }
                        SKL[k][l] += temp * (float)Nv[i][degreeV - l];
                    }
                }
            }

            return SKL;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetRationalDerivatives(int, Vector2)"/>
        public static Vector2[][] ComputeRationalSurfaceDerivatives(
            in NurbsSurface<Vector2> surface,
            int derivative,
            Vector2 uv
        )
        {
            var degreeU = surface.DegreeU;
            var degreeV = surface.DegreeV;
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;
            var controlPoints = surface.ControlPoints;

            Validate.Argument(
                derivative > 0,
                nameof(derivative),
                "Derivative must be greater than zero."
            );
            Validate.Range(uv.X, knotsU[0], knotsU[knotsU.Count - 1], nameof(uv.X));
            Validate.Range(uv.Y, knotsV[0], knotsV[knotsV.Count - 1], nameof(uv.Y));

            // Aders = weighted position derivatives (reuse existing ComputeDerivatives)
            var Aders = ComputeDerivatives(in surface, derivative, uv);

            // Compute wders (weight derivatives) separately
            int uSpanIndex = Polynomials.GetKnotSpanIndex(degreeU, knotsU, uv.X);
            double[][] Nu = Polynomials.BasisFunctionsDerivatives(
                uSpanIndex,
                degreeU,
                derivative,
                knotsU,
                uv.X
            );

            int vSpanIndex = Polynomials.GetKnotSpanIndex(degreeV, knotsV, uv.Y);
            double[][] Nv = Polynomials.BasisFunctionsDerivatives(
                vSpanIndex,
                degreeV,
                derivative,
                knotsV,
                uv.Y
            );

            int du = Math.Min(derivative, degreeU);
            int dv = Math.Min(derivative, degreeV);

            var wders = new double[derivative + 1][];
            for (int i = 0; i <= derivative; i++)
            {
                wders[i] = new double[derivative + 1];
            }

            var tempW = new double[degreeV + 1];

            for (int k = 0; k <= du; k++)
            {
                for (int s = 0; s <= degreeV; s++)
                {
                    tempW[s] = 0.0;
                    for (int r = 0; r <= degreeU; r++)
                    {
                        var cp = controlPoints[uSpanIndex - degreeU + r][vSpanIndex - degreeV + s];
                        tempW[s] += cp.Weight * Nu[k][r];
                    }
                }

                int dd = Math.Min(derivative - k, dv);
                for (int l = 0; l <= dd; l++)
                {
                    for (int s = 0; s <= degreeV; s++)
                    {
                        wders[k][l] += tempW[s] * Nv[l][s];
                    }
                }
            }

            // Apply quotient rule for rational derivatives
            var derivatives = new Vector2[derivative + 1][];
            for (int i = 0; i <= derivative; i++)
            {
                derivatives[i] = new Vector2[derivative + 1];
            }

            for (int k = 0; k <= derivative; k++)
            {
                for (int l = 0; l <= derivative - k; l++)
                {
                    var v = Aders[k][l];

                    for (int j = 1; j <= l; j++)
                    {
                        v -=
                            (float)(MathUtils.Binomial(l, j) * wders[0][j]) * derivatives[k][l - j];
                    }

                    for (int i = 1; i <= k; i++)
                    {
                        v -=
                            (float)(MathUtils.Binomial(k, i) * wders[i][0]) * derivatives[k - i][l];

                        var v2 = Vector2.Zero;
                        for (int j = 1; j <= l; j++)
                        {
                            v2 +=
                                (float)(MathUtils.Binomial(l, j) * wders[i][j])
                                * derivatives[k - i][l - j];
                        }
                        v -= (float)MathUtils.Binomial(k, i) * v2;
                    }

                    derivatives[k][l] = MathUtils.IsZero(wders[0][0])
                        ? Vector2.Zero
                        : v / (float)wders[0][0];
                }
            }

            return derivatives;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetRationalFirstOrderDerivatives(Vector2, out T, out T, out T)"/>
        public static void ComputeRationalSurfaceFirstOrderDerivatives(
            in NurbsSurface<Vector2> surface,
            Vector2 uv,
            out Vector2 S,
            out Vector2 Su,
            out Vector2 Sv
        )
        {
            var degreeU = surface.DegreeU;
            var degreeV = surface.DegreeV;
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;
            var controlPoints = surface.ControlPoints;

            Validate.Range(uv.X, knotsU[0], knotsU[knotsU.Count - 1], nameof(uv.X));
            Validate.Range(uv.Y, knotsV[0], knotsV[knotsV.Count - 1], nameof(uv.Y));

            int uSpanIndex = Polynomials.GetKnotSpanIndex(degreeU, knotsU, uv.X);
            double[][] Nu = Polynomials.BasisFunctionsFirstOrderDerivative(
                uSpanIndex,
                degreeU,
                knotsU,
                uv.X
            );

            int vSpanIndex = Polynomials.GetKnotSpanIndex(degreeV, knotsV, uv.Y);
            double[][] Nv = Polynomials.BasisFunctionsFirstOrderDerivative(
                vSpanIndex,
                degreeV,
                knotsV,
                uv.Y
            );

            int du = Math.Min(1, degreeU);
            int dv = Math.Min(1, degreeV);

            // Track Aders (weighted position) and wders (weight) for 2x2
            var Aders = new Vector2[2, 2];
            var wders = new double[2, 2];
            var tempA = new Vector2[degreeV + 1];
            var tempW = new double[degreeV + 1];

            for (int k = 0; k <= du; k++)
            {
                for (int s = 0; s <= degreeV; s++)
                {
                    tempA[s] = Vector2.Zero;
                    tempW[s] = 0.0;
                    for (int r = 0; r <= degreeU; r++)
                    {
                        var cp = controlPoints[uSpanIndex - degreeU + r][vSpanIndex - degreeV + s];
                        double wBasis = cp.Weight * Nu[k][r];
                        tempA[s] += cp.Value * (float)wBasis;
                        tempW[s] += wBasis;
                    }
                }

                int dd = Math.Min(1 - k, dv);
                for (int l = 0; l <= dd; l++)
                {
                    for (int s = 0; s <= degreeV; s++)
                    {
                        Aders[k, l] += tempA[s] * (float)Nv[l][s];
                        wders[k, l] += tempW[s] * Nv[l][s];
                    }
                }
            }

            // Apply quotient rule
            double invW = MathUtils.IsZero(wders[0, 0]) ? 0.0 : 1.0 / wders[0, 0];
            S = Aders[0, 0] * (float)invW;
            Sv = (Aders[0, 1] - (float)wders[0, 1] * S) * (float)invW;
            Su = (Aders[1, 0] - (float)wders[1, 0] * S) * (float)invW;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetNormal(Vector2)"/>
        public static Vector2 ComputeNormal(in NurbsSurface<Vector2> surface, Vector2 uv)
        {
            return Vector2.Zero;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetCurvature(SurfaceCurvature, Vector2)"/>
        public static double ComputeCurvature(
            in NurbsSurface<Vector2> surface,
            SurfaceCurvature curvature,
            Vector2 uv
        )
        {
            return 0.0;
        }

        /// <inheritdoc cref="NurbsSurface{T}.SwapDim()"/>
        public static NurbsSurface<Vector2> SwapDim(in NurbsSurface<Vector2> surface)
        {
            var controlPoints = surface.ControlPoints;
            int countU = controlPoints.Count;
            int countV = controlPoints[0].Count;

            // Transpose control points
            var transposedControlPoints = new ControlPoint<Vector2>[countV][];
            for (int j = 0; j < countV; j++)
            {
                transposedControlPoints[j] = new ControlPoint<Vector2>[countU];
                for (int i = 0; i < countU; i++)
                {
                    transposedControlPoints[j][i] = controlPoints[i][j];
                }
            }

            return new NurbsSurface<Vector2>(
                surface.DegreeV, // New DegreeU = old DegreeV
                surface.DegreeU, // New DegreeV = old DegreeU
                transposedControlPoints,
                surface.KnotsV, // New KnotsU = old KnotsV
                surface.KnotsU // New KnotsV = old KnotsU
            );
        }

        /// <inheritdoc cref="NurbsSurface{T}.Reverse(SurfaceDirection)"/>
        public static NurbsSurface<Vector2> Reverse(
            in NurbsSurface<Vector2> surface,
            SurfaceDirection direction
        )
        {
            var controlPoints = surface.ControlPoints;
            int countU = controlPoints.Count;
            int countV = controlPoints[0].Count;

            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;

            // Reverse knot vectors as needed
            IReadOnlyList<double> newKnotsU =
                (direction == SurfaceDirection.All || direction == SurfaceDirection.UDirection)
                    ? KnotsUtils.ReverseKnots(knotsU)
                    : knotsU;

            IReadOnlyList<double> newKnotsV =
                (direction == SurfaceDirection.All || direction == SurfaceDirection.VDirection)
                    ? KnotsUtils.ReverseKnots(knotsV)
                    : knotsV;

            // Reverse control points
            var newControlPoints = new ControlPoint<Vector2>[countU][];

            if (direction == SurfaceDirection.UDirection)
            {
                // Reverse each column (swap rows for each column index)
                for (int i = 0; i < countU; i++)
                {
                    newControlPoints[i] = new ControlPoint<Vector2>[countV];
                    for (int j = 0; j < countV; j++)
                    {
                        newControlPoints[i][j] = controlPoints[countU - 1 - i][j];
                    }
                }
            }
            else if (direction == SurfaceDirection.VDirection)
            {
                // Reverse each row
                for (int i = 0; i < countU; i++)
                {
                    newControlPoints[i] = new ControlPoint<Vector2>[countV];
                    for (int j = 0; j < countV; j++)
                    {
                        newControlPoints[i][j] = controlPoints[i][countV - 1 - j];
                    }
                }
            }
            else // All
            {
                // Reverse both directions
                for (int i = 0; i < countU; i++)
                {
                    newControlPoints[i] = new ControlPoint<Vector2>[countV];
                    for (int j = 0; j < countV; j++)
                    {
                        newControlPoints[i][j] = controlPoints[countU - 1 - i][countV - 1 - j];
                    }
                }
            }

            return new NurbsSurface<Vector2>(
                surface.DegreeU,
                surface.DegreeV,
                newControlPoints,
                newKnotsU,
                newKnotsV
            );
        }

        /// <inheritdoc cref="NurbsSurface{T}.InsertKnot(double, int, SurfaceDirection, out NurbsSurface{T})"/>
        public static int InsertKnot(
            in NurbsSurface<Vector2> surface,
            double insertKnot,
            int times,
            SurfaceDirection direction,
            out NurbsSurface<Vector2> result
        )
        {
            Validate.Argument(
                direction == SurfaceDirection.UDirection
                    || direction == SurfaceDirection.VDirection,
                nameof(direction),
                "Direction must be UDirection or VDirection."
            );
            Validate.Argument(times > 0, nameof(times), "Times must be greater than zero.");

            bool isUDirection = direction == SurfaceDirection.UDirection;
            int degree = isUDirection ? surface.DegreeU : surface.DegreeV;
            var knotVector = isUDirection ? surface.KnotsU : surface.KnotsV;
            var controlPoints = surface.ControlPoints;

            int knotSpanIndex = Polynomials.GetKnotSpanIndex(degree, knotVector, insertKnot);
            int multiplicity = Polynomials.GetKnotMultiplicity(knotVector, insertKnot);

            if (multiplicity >= degree)
            {
                result = surface;
                return 0;
            }

            if (times + multiplicity > degree)
            {
                times = degree - multiplicity;
            }

            if (times <= 0)
            {
                result = surface;
                return 0;
            }

            // Build new knot vector
            var insertedKnotVector = new double[knotVector.Count + times];
            for (int i = 0; i <= knotSpanIndex; i++)
            {
                insertedKnotVector[i] = knotVector[i];
            }
            for (int i = 1; i <= times; i++)
            {
                insertedKnotVector[knotSpanIndex + i] = insertKnot;
            }
            for (int i = knotSpanIndex + 1; i < knotVector.Count; i++)
            {
                insertedKnotVector[i + times] = knotVector[i];
            }

            // Compute alpha coefficients
            var alpha = new double[degree - multiplicity][];
            for (int i = 0; i < degree - multiplicity; i++)
            {
                alpha[i] = new double[times + 1];
            }

            for (int j = 1; j <= times; j++)
            {
                int L = knotSpanIndex - degree + j;
                for (int i = 0; i <= degree - j - multiplicity; i++)
                {
                    alpha[i][j] =
                        (insertKnot - knotVector[L + i])
                        / (knotVector[i + knotSpanIndex + 1] - knotVector[L + i]);
                }
            }

            var temp = new ControlPoint<Vector2>[degree + 1];
            int rows = controlPoints.Count;
            int columns = controlPoints[0].Count;

            if (isUDirection)
            {
                var updatedControlPoints = new ControlPoint<Vector2>[rows + times][];
                for (int i = 0; i < rows + times; i++)
                {
                    updatedControlPoints[i] = new ControlPoint<Vector2>[columns];
                }

                for (int col = 0; col < columns; col++)
                {
                    // Copy unaffected control points
                    for (int i = 0; i <= knotSpanIndex - degree; i++)
                    {
                        updatedControlPoints[i][col] = controlPoints[i][col];
                    }
                    for (int i = knotSpanIndex - multiplicity; i < rows; i++)
                    {
                        updatedControlPoints[i + times][col] = controlPoints[i][col];
                    }

                    // Load affected control points into temp
                    for (int i = 0; i < degree - multiplicity + 1; i++)
                    {
                        temp[i] = controlPoints[knotSpanIndex - degree + i][col];
                    }

                    // Insert knot
                    int L = 0;
                    for (int j = 1; j <= times; j++)
                    {
                        L = knotSpanIndex - degree + j;
                        for (int i = 0; i <= degree - j - multiplicity; i++)
                        {
                            double a = alpha[i][j];
                            temp[i] = ControlPointsHelper.BlendControlPoints(
                                temp[i],
                                temp[i + 1],
                                a
                            );
                        }
                        updatedControlPoints[L][col] = temp[0];
                        updatedControlPoints[knotSpanIndex + times - j - multiplicity][col] = temp[
                            degree - j - multiplicity
                        ];
                    }

                    for (int i = L + 1; i < knotSpanIndex - multiplicity; i++)
                    {
                        updatedControlPoints[i][col] = temp[i - L];
                    }
                }

                result = new NurbsSurface<Vector2>(
                    surface.DegreeU,
                    surface.DegreeV,
                    updatedControlPoints,
                    insertedKnotVector,
                    surface.KnotsV
                );
            }
            else
            {
                var updatedControlPoints = new ControlPoint<Vector2>[rows][];
                for (int i = 0; i < rows; i++)
                {
                    updatedControlPoints[i] = new ControlPoint<Vector2>[columns + times];
                }

                for (int row = 0; row < rows; row++)
                {
                    // Copy unaffected control points
                    for (int i = 0; i <= knotSpanIndex - degree; i++)
                    {
                        updatedControlPoints[row][i] = controlPoints[row][i];
                    }
                    for (int i = knotSpanIndex - multiplicity; i < columns; i++)
                    {
                        updatedControlPoints[row][i + times] = controlPoints[row][i];
                    }

                    // Load affected control points into temp
                    for (int i = 0; i < degree - multiplicity + 1; i++)
                    {
                        temp[i] = controlPoints[row][knotSpanIndex - degree + i];
                    }

                    // Insert knot
                    int L = 0;
                    for (int j = 1; j <= times; j++)
                    {
                        L = knotSpanIndex - degree + j;
                        for (int i = 0; i <= degree - j - multiplicity; i++)
                        {
                            double a = alpha[i][j];
                            temp[i] = ControlPointsHelper.BlendControlPoints(
                                temp[i],
                                temp[i + 1],
                                a
                            );
                        }
                        updatedControlPoints[row][L] = temp[0];
                        updatedControlPoints[row][knotSpanIndex + times - j - multiplicity] = temp[
                            degree - j - multiplicity
                        ];
                    }

                    for (int i = L + 1; i < knotSpanIndex - multiplicity; i++)
                    {
                        updatedControlPoints[row][i] = temp[i - L];
                    }
                }

                result = new NurbsSurface<Vector2>(
                    surface.DegreeU,
                    surface.DegreeV,
                    updatedControlPoints,
                    surface.KnotsU,
                    insertedKnotVector
                );
            }

            return times;
        }

        /// <inheritdoc cref="NurbsSurface{T}.RefineKnotVector(IReadOnlyList{double}, SurfaceDirection)"/>
        public static NurbsSurface<Vector2> RefineKnotVector(
            in NurbsSurface<Vector2> surface,
            IReadOnlyList<double> insertKnotElements,
            SurfaceDirection direction
        )
        {
            Validate.Argument(
                direction == SurfaceDirection.UDirection
                    || direction == SurfaceDirection.VDirection,
                nameof(direction),
                "Direction must be UDirection or VDirection."
            );
            Validate.Argument(
                insertKnotElements != null && insertKnotElements.Count > 0,
                nameof(insertKnotElements),
                "insertKnotElements size must be greater than zero."
            );

            var controlPoints = surface.ControlPoints;
            bool isUDirection = direction == SurfaceDirection.UDirection;

            if (isUDirection)
            {
                // Transpose, refine each row, transpose back
                var transposed = MathUtils.Transpose(controlPoints);
                var tempControlPoints = new List<ControlPoint<Vector2>[]>();
                IReadOnlyList<double> newKnotsU = null;

                for (int i = 0; i < transposed.Count; i++)
                {
                    var curve = new NurbsCurve<Vector2>(
                        surface.DegreeU,
                        transposed[i],
                        surface.KnotsU
                    );
                    var refined = NurbsCurve2DHelper.RefineKnotVector(curve, insertKnotElements);

                    var cpArray = new ControlPoint<Vector2>[refined.ControlPoints.Count];
                    for (int j = 0; j < refined.ControlPoints.Count; j++)
                    {
                        cpArray[j] = refined.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsU = refined.Knots;
                }

                // Convert to array and transpose back
                var tempArray = tempControlPoints.ToArray();
                var updatedControlPoints = MathUtils.Transpose(tempArray);

                return new NurbsSurface<Vector2>(
                    surface.DegreeU,
                    surface.DegreeV,
                    updatedControlPoints,
                    newKnotsU,
                    surface.KnotsV
                );
            }
            else
            {
                // Refine each row directly
                var tempControlPoints = new List<ControlPoint<Vector2>[]>();
                IReadOnlyList<double> newKnotsV = null;

                for (int i = 0; i < controlPoints.Count; i++)
                {
                    var rowCPs = new ControlPoint<Vector2>[controlPoints[i].Count];
                    for (int j = 0; j < controlPoints[i].Count; j++)
                    {
                        rowCPs[j] = controlPoints[i][j];
                    }

                    var curve = new NurbsCurve<Vector2>(surface.DegreeV, rowCPs, surface.KnotsV);
                    var refined = NurbsCurve2DHelper.RefineKnotVector(curve, insertKnotElements);

                    var cpArray = new ControlPoint<Vector2>[refined.ControlPoints.Count];
                    for (int j = 0; j < refined.ControlPoints.Count; j++)
                    {
                        cpArray[j] = refined.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsV = refined.Knots;
                }

                return new NurbsSurface<Vector2>(
                    surface.DegreeU,
                    surface.DegreeV,
                    tempControlPoints.ToArray(),
                    surface.KnotsU,
                    newKnotsV
                );
            }
        }

        /// <inheritdoc cref="NurbsSurface{T}.DecomposeToBeziers()"/>
        public static IReadOnlyList<BezierSurface<Vector2>> DecomposeToBeziers(
            in NurbsSurface<Vector2> surface
        )
        {
            int degreeU = surface.DegreeU;
            int degreeV = surface.DegreeV;
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;
            var controlPoints = surface.ControlPoints;

            int rows = controlPoints.Count;
            int columns = controlPoints[0].Count;

            // Phase 1: Decompose in U direction
            int tempPatchCount = rows - degreeU;
            var tempPatches = new ControlPoint<Vector2>[tempPatchCount][][];
            for (int i = 0; i < tempPatchCount; i++)
            {
                tempPatches[i] = new ControlPoint<Vector2>[degreeU + 1][];
                for (int j = 0; j <= degreeU; j++)
                {
                    tempPatches[i][j] = new ControlPoint<Vector2>[columns];
                }
            }

            int m = rows - 1 + degreeU + 1;
            int a = degreeU;
            int b = degreeU + 1;
            int nb = 0;

            // Copy first patch control points
            for (int i = 0; i <= degreeU; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    tempPatches[nb][i][j] = controlPoints[i][j];
                }
            }

            var alphaVector = new double[Math.Max(degreeU, degreeV) + 1];

            while (b < m)
            {
                int ii = b;
                while (b < m && MathUtils.IsAlmostEqualTo(knotsU[b + 1], knotsU[b]))
                {
                    b++;
                }
                int multi = b - ii + 1;

                if (multi < degreeU)
                {
                    double numerator = knotsU[b] - knotsU[a];

                    for (int j = degreeU; j > multi; j--)
                    {
                        alphaVector[j - multi - 1] = numerator / (knotsU[a + j] - knotsU[a]);
                    }

                    int r = degreeU - multi;
                    for (int j = 1; j <= r; j++)
                    {
                        int save = r - j;
                        int s = multi + j;

                        for (int k = degreeU; k >= s; k--)
                        {
                            double alpha = alphaVector[k - s];
                            for (int col = 0; col < columns; col++)
                            {
                                tempPatches[nb][k][col] = ControlPointsHelper.BlendControlPoints(
                                    tempPatches[nb][k - 1][col],
                                    tempPatches[nb][k][col],
                                    alpha
                                );
                            }
                        }

                        if (b < m && nb + 1 < tempPatchCount)
                        {
                            for (int col = 0; col < columns; col++)
                            {
                                tempPatches[nb + 1][save][col] = tempPatches[nb][degreeU][col];
                            }
                        }
                    }
                }

                nb++;
                if (b < m && nb < tempPatchCount)
                {
                    for (int i = degreeU - multi; i <= degreeU; i++)
                    {
                        for (int col = 0; col < columns; col++)
                        {
                            tempPatches[nb][i][col] = controlPoints[b - degreeU + i][col];
                        }
                    }
                    a = b;
                    b++;
                }
            }

            int tempSize = nb;

            // Phase 2: Decompose each temp patch in V direction
            int finalPatchCount = tempSize * (columns - degreeV);
            var bezierPatches = new List<ControlPoint<Vector2>[][]>(finalPatchCount);

            for (int i = 0; i < finalPatchCount; i++)
            {
                var patch = new ControlPoint<Vector2>[degreeU + 1][];
                for (int j = 0; j <= degreeU; j++)
                {
                    patch[j] = new ControlPoint<Vector2>[degreeV + 1];
                }
                bezierPatches.Add(patch);
            }

            nb = 0;
            for (int np = 0; np < tempSize; np++)
            {
                // Copy initial patch
                for (int i = 0; i <= degreeU; i++)
                {
                    for (int j = 0; j <= degreeV; j++)
                    {
                        bezierPatches[nb][i][j] = tempPatches[np][i][j];
                    }
                }

                m = columns + degreeV;
                a = degreeV;
                b = degreeV + 1;

                while (b < m)
                {
                    int ii = b;
                    while (b < m && MathUtils.IsAlmostEqualTo(knotsV[b + 1], knotsV[b]))
                    {
                        b++;
                    }
                    int multi = b - ii + 1;

                    if (multi < degreeV)
                    {
                        double numerator = knotsV[b] - knotsV[a];

                        for (int j = degreeV; j > multi; j--)
                        {
                            alphaVector[j - multi - 1] = numerator / (knotsV[a + j] - knotsV[a]);
                        }

                        int r = degreeV - multi;
                        for (int j = 1; j <= r; j++)
                        {
                            int save = r - j;
                            int s = multi + j;

                            for (int k = degreeV; k >= s; k--)
                            {
                                double alpha = alphaVector[k - s];
                                for (int row = 0; row <= degreeU; row++)
                                {
                                    bezierPatches[nb][row][k] =
                                        ControlPointsHelper.BlendControlPoints(
                                            bezierPatches[nb][row][k - 1],
                                            bezierPatches[nb][row][k],
                                            alpha
                                        );
                                }
                            }

                            if (b < m && nb + 1 < bezierPatches.Count)
                            {
                                for (int row = 0; row <= degreeU; row++)
                                {
                                    bezierPatches[nb + 1][row][save] = bezierPatches[nb][row][
                                        degreeV
                                    ];
                                }
                            }
                        }
                    }

                    nb++;
                    if (b < m && nb < bezierPatches.Count)
                    {
                        for (int i = degreeV - multi; i <= degreeV; i++)
                        {
                            for (int row = 0; row <= degreeU; row++)
                            {
                                bezierPatches[nb][row][i] = tempPatches[np][row][b - degreeV + i];
                            }
                        }
                        a = b;
                        b++;
                    }
                }
            }

            // Trim excess patches and convert to BezierSurface
            var result = new List<BezierSurface<Vector2>>(nb);
            for (int i = 0; i < nb; i++)
            {
                result.Add(new BezierSurface<Vector2>(degreeU, degreeV, bezierPatches[i]));
            }

            return result;
        }

        /// <inheritdoc cref="NurbsSurface{T}.TryRemoveKnot(double, int, SurfaceDirection, out NurbsSurface{T})"/>`
        public static bool RemoveKnot(
            in NurbsSurface<Vector2> surface,
            double removeKnot,
            int times,
            SurfaceDirection direction,
            out NurbsSurface<Vector2> result
        )
        {
            Validate.Argument(
                direction == SurfaceDirection.UDirection
                    || direction == SurfaceDirection.VDirection,
                nameof(direction),
                "Direction must be UDirection or VDirection."
            );
            Validate.Argument(times > 0, nameof(times), "Times must be greater than zero.");

            var controlPoints = surface.ControlPoints;
            bool isUDirection = direction == SurfaceDirection.UDirection;

            if (isUDirection)
            {
                Validate.Range(
                    removeKnot,
                    surface.KnotsU[0],
                    surface.KnotsU[surface.KnotsU.Count - 1],
                    nameof(removeKnot)
                );

                // Transpose, remove knot from each row, transpose back
                var transposed = MathUtils.Transpose(controlPoints);
                var tempControlPoints = new List<ControlPoint<Vector2>[]>();
                IReadOnlyList<double> newKnotsU = null;

                for (int i = 0; i < transposed.Count; i++)
                {
                    var curve = new NurbsCurve<Vector2>(
                        surface.DegreeU,
                        transposed[i],
                        surface.KnotsU
                    );
                    if (!NurbsCurve2DHelper.RemoveKnot(curve, removeKnot, times, out var removed))
                    {
                        result = surface;
                        return false;
                    }

                    var cpArray = new ControlPoint<Vector2>[removed.ControlPoints.Count];
                    for (int j = 0; j < removed.ControlPoints.Count; j++)
                    {
                        cpArray[j] = removed.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsU = removed.Knots;
                }

                var tempArray = tempControlPoints.ToArray();
                var updatedControlPoints = MathUtils.Transpose(tempArray);

                result = new NurbsSurface<Vector2>(
                    surface.DegreeU,
                    surface.DegreeV,
                    updatedControlPoints,
                    newKnotsU,
                    surface.KnotsV
                );
            }
            else
            {
                Validate.Range(
                    removeKnot,
                    surface.KnotsV[0],
                    surface.KnotsV[surface.KnotsV.Count - 1],
                    nameof(removeKnot)
                );

                var tempControlPoints = new List<ControlPoint<Vector2>[]>();
                IReadOnlyList<double> newKnotsV = null;

                for (int i = 0; i < controlPoints.Count; i++)
                {
                    var rowCPs = new ControlPoint<Vector2>[controlPoints[i].Count];
                    for (int j = 0; j < controlPoints[i].Count; j++)
                    {
                        rowCPs[j] = controlPoints[i][j];
                    }

                    var curve = new NurbsCurve<Vector2>(surface.DegreeV, rowCPs, surface.KnotsV);
                    if (!NurbsCurve2DHelper.RemoveKnot(curve, removeKnot, times, out var removed))
                    {
                        result = surface;
                        return false;
                    }

                    var cpArray = new ControlPoint<Vector2>[removed.ControlPoints.Count];
                    for (int j = 0; j < removed.ControlPoints.Count; j++)
                    {
                        cpArray[j] = removed.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsV = removed.Knots;
                }

                result = new NurbsSurface<Vector2>(
                    surface.DegreeU,
                    surface.DegreeV,
                    tempControlPoints.ToArray(),
                    surface.KnotsU,
                    newKnotsV
                );
            }

            return true;
        }

        public static NurbsSurface<Vector2> ElevateDegree(
            in NurbsSurface<Vector2> surface,
            int times,
            SurfaceDirection direction
        )
        {
            Validate.Argument(
                direction == SurfaceDirection.UDirection
                    || direction == SurfaceDirection.VDirection,
                nameof(direction),
                "Direction must be UDirection or VDirection."
            );
            Validate.Argument(times > 0, nameof(times), "Times must be greater than zero.");

            var controlPoints = surface.ControlPoints;
            bool isUDirection = direction == SurfaceDirection.UDirection;

            if (isUDirection)
            {
                // Transpose, elevate degree for each row, transpose back
                var transposed = MathUtils.Transpose(controlPoints);
                var tempControlPoints = new List<ControlPoint<Vector2>[]>();
                IReadOnlyList<double> newKnotsU = null;

                for (int i = 0; i < transposed.Count; i++)
                {
                    var curve = new NurbsCurve<Vector2>(
                        surface.DegreeU,
                        transposed[i],
                        surface.KnotsU
                    );
                    var elevated = NurbsCurve2DHelper.ElevateDegree(curve, times);

                    var cpArray = new ControlPoint<Vector2>[elevated.ControlPoints.Count];
                    for (int j = 0; j < elevated.ControlPoints.Count; j++)
                    {
                        cpArray[j] = elevated.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsU = elevated.Knots;
                }

                var tempArray = tempControlPoints.ToArray();
                var updatedControlPoints = MathUtils.Transpose(tempArray);

                return new NurbsSurface<Vector2>(
                    surface.DegreeU + times,
                    surface.DegreeV,
                    updatedControlPoints,
                    newKnotsU,
                    surface.KnotsV
                );
            }
            else
            {
                // Elevate degree for each row directly
                var tempControlPoints = new List<ControlPoint<Vector2>[]>();
                IReadOnlyList<double> newKnotsV = null;

                for (int i = 0; i < controlPoints.Count; i++)
                {
                    var rowCPs = new ControlPoint<Vector2>[controlPoints[i].Count];
                    for (int j = 0; j < controlPoints[i].Count; j++)
                    {
                        rowCPs[j] = controlPoints[i][j];
                    }

                    var curve = new NurbsCurve<Vector2>(surface.DegreeV, rowCPs, surface.KnotsV);
                    var elevated = NurbsCurve2DHelper.ElevateDegree(curve, times);

                    var cpArray = new ControlPoint<Vector2>[elevated.ControlPoints.Count];
                    for (int j = 0; j < elevated.ControlPoints.Count; j++)
                    {
                        cpArray[j] = elevated.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsV = elevated.Knots;
                }

                return new NurbsSurface<Vector2>(
                    surface.DegreeU,
                    surface.DegreeV + times,
                    tempControlPoints.ToArray(),
                    surface.KnotsU,
                    newKnotsV
                );
            }
        }

        public static bool ReduceDegree(
            in NurbsSurface<Vector2> surface,
            SurfaceDirection direction,
            out NurbsSurface<Vector2> result
        )
        {
            Validate.Argument(
                direction == SurfaceDirection.UDirection
                    || direction == SurfaceDirection.VDirection,
                nameof(direction),
                "Direction must be UDirection or VDirection."
            );

            var controlPoints = surface.ControlPoints;
            bool isUDirection = direction == SurfaceDirection.UDirection;

            if (isUDirection)
            {
                if (surface.DegreeU <= 1)
                {
                    result = surface;
                    return false;
                }

                // Transpose, reduce degree for each row, transpose back
                var transposed = MathUtils.Transpose(controlPoints);
                var tempControlPoints = new List<ControlPoint<Vector2>[]>();
                IReadOnlyList<double> newKnotsU = null;

                for (int i = 0; i < transposed.Length; i++)
                {
                    var curve = new NurbsCurve<Vector2>(
                        surface.DegreeU,
                        transposed[i],
                        surface.KnotsU
                    );
                    if (!NurbsCurve2DHelper.ReduceDegree(curve, out var reduced))
                    {
                        result = surface;
                        return false;
                    }

                    var cpArray = new ControlPoint<Vector2>[reduced.ControlPoints.Count];
                    for (int j = 0; j < reduced.ControlPoints.Count; j++)
                    {
                        cpArray[j] = reduced.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsU = reduced.Knots;
                }

                var tempArray = tempControlPoints.ToArray();
                var updatedControlPoints = MathUtils.Transpose(tempArray);

                result = new NurbsSurface<Vector2>(
                    surface.DegreeU - 1,
                    surface.DegreeV,
                    updatedControlPoints,
                    newKnotsU,
                    surface.KnotsV
                );
            }
            else
            {
                if (surface.DegreeV <= 1)
                {
                    result = surface;
                    return false;
                }

                // Reduce degree for each row directly
                var tempControlPoints = new List<ControlPoint<Vector2>[]>();
                IReadOnlyList<double> newKnotsV = null;

                for (int i = 0; i < controlPoints.Count; i++)
                {
                    var rowCPs = new ControlPoint<Vector2>[controlPoints[i].Count];
                    for (int j = 0; j < controlPoints[i].Count; j++)
                    {
                        rowCPs[j] = controlPoints[i][j];
                    }

                    var curve = new NurbsCurve<Vector2>(surface.DegreeV, rowCPs, surface.KnotsV);
                    if (!NurbsCurve2DHelper.ReduceDegree(curve, out var reduced))
                    {
                        result = surface;
                        return false;
                    }

                    var cpArray = new ControlPoint<Vector2>[reduced.ControlPoints.Count];
                    for (int j = 0; j < reduced.ControlPoints.Count; j++)
                    {
                        cpArray[j] = reduced.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsV = reduced.Knots;
                }

                result = new NurbsSurface<Vector2>(
                    surface.DegreeU,
                    surface.DegreeV - 1,
                    tempControlPoints.ToArray(),
                    surface.KnotsU,
                    newKnotsV
                );
            }

            return true;
        }
    }
}
