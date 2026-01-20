using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    public static class NurbsSurface3DHelper
    {
        /// <inheritdoc cref="NurbsSurface{T}.GetPointOnSurface(Vector2)"/>
        public static Vector3 GetPointOnSurface(in NurbsSurface<Vector3> surface, Vector2 uv)
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
            var result = new ControlPoint<Vector3>(Vector3.Zero, 0.0);

            for (int l = 0; l <= degreeV; l++)
            {
                var temp = new ControlPoint<Vector3>(Vector3.Zero, 0.0);
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

            return result.Weight != 0 ? result.Value / (float)result.Weight : Vector3.Zero;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetDerivatives(int, Vector2)"/>
        public static Vector3[][] ComputeDerivatives(
            in NurbsSurface<Vector3> surface,
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
            var derivatives = new Vector3[derivative + 1][];
            for (int i = 0; i <= derivative; i++)
            {
                derivatives[i] = new Vector3[derivative + 1];
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

            var temp = new Vector3[degreeV + 1];

            for (int k = 0; k <= du; k++)
            {
                for (int s = 0; s <= degreeV; s++)
                {
                    temp[s] = Vector3.Zero;
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
        public static Vector3[][] ComputeFirstOrderDerivatives(
            in NurbsSurface<Vector3> surface,
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

            var derivatives = new Vector3[2][] { new Vector3[2], new Vector3[2] };

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

            var temp = new Vector3[degreeV + 1];

            for (int k = 0; k <= du; k++)
            {
                for (int s = 0; s <= degreeV; s++)
                {
                    temp[s] = Vector3.Zero;
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
        public static Vector3[][][][] ComputeControlPointsOfDerivatives(
            in NurbsSurface<Vector3> surface,
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
            var PKL = new Vector3[derivative + 1][][][];
            for (int k = 0; k <= derivative; k++)
            {
                PKL[k] = new Vector3[derivative + 1][][];
                for (int l = 0; l <= derivative; l++)
                {
                    PKL[k][l] = new Vector3[rangeU + 1][];
                    for (int i = 0; i <= rangeU; i++)
                    {
                        PKL[k][l][i] = new Vector3[rangeV + 1];
                    }
                }
            }

            // Compute derivatives in U direction for each V column
            for (int j = minSpanIndexV; j <= maxSpanIndexV; j++)
            {
                var points = new ControlPoint<Vector3>[controlPoints.Count];
                for (int i = 0; i < controlPoints.Count; i++)
                {
                    points[i] = controlPoints[i][j];
                }

                var tempCurve = new NurbsCurve<Vector3>(degreeU, points, knotsU);
                var temp = NurbsCurve3DHelper.ComputeControlPointsOfDerivatives(
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

                    var points = new ControlPoint<Vector3>[rangeV + 1];
                    for (int j = 0; j <= rangeV; j++)
                    {
                        points[j] = new ControlPoint<Vector3>(PKL[k][0][i][j], 1.0);
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

                    var tempCurve = new NurbsCurve<Vector3>(degreeV, points, tempKnots);
                    var temp = NurbsCurve3DHelper.ComputeControlPointsOfDerivatives(
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
        public static Vector3[][] ComputeDerivativesByAllBasisFunctions(
            in NurbsSurface<Vector3> surface,
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

            var SKL = new Vector3[derivative + 1][];
            for (int i = 0; i <= derivative; i++)
            {
                SKL[i] = new Vector3[derivative + 1];
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
                    SKL[k][l] = Vector3.Zero;
                    for (int i = 0; i <= degreeV - l; i++)
                    {
                        var temp = Vector3.Zero;
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
        public static Vector3[][] ComputeRationalSurfaceDerivatives(
            in NurbsSurface<Vector3> surface,
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
            var derivatives = new Vector3[derivative + 1][];
            for (int i = 0; i <= derivative; i++)
            {
                derivatives[i] = new Vector3[derivative + 1];
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

                        var v2 = Vector3.Zero;
                        for (int j = 1; j <= l; j++)
                        {
                            v2 +=
                                (float)(MathUtils.Binomial(l, j) * wders[i][j])
                                * derivatives[k - i][l - j];
                        }
                        v -= (float)MathUtils.Binomial(k, i) * v2;
                    }

                    derivatives[k][l] = MathUtils.IsZero(wders[0][0])
                        ? Vector3.Zero
                        : v / (float)wders[0][0];
                }
            }

            return derivatives;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetRationalFirstOrderDerivatives(Vector2, out T, out T, out T)"/>
        public static void ComputeRationalSurfaceFirstOrderDerivatives(
            in NurbsSurface<Vector3> surface,
            Vector2 uv,
            out Vector3 S,
            out Vector3 Su,
            out Vector3 Sv
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
            var Aders = new Vector3[2, 2];
            var wders = new double[2, 2];
            var tempA = new Vector3[degreeV + 1];
            var tempW = new double[degreeV + 1];

            for (int k = 0; k <= du; k++)
            {
                for (int s = 0; s <= degreeV; s++)
                {
                    tempA[s] = Vector3.Zero;
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
        public static Vector3 ComputeNormal(in NurbsSurface<Vector3> surface, Vector2 uv)
        {
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;

            Validate.Range(uv.X, knotsU[0], knotsU[knotsU.Count - 1], nameof(uv.X));
            Validate.Range(uv.Y, knotsV[0], knotsV[knotsV.Count - 1], nameof(uv.Y));

            var derivatives = ComputeRationalSurfaceDerivatives(in surface, 1, uv);
            var Su = derivatives[1][0];
            var Sv = derivatives[0][1];

            // Normalize Su, then cross with Sv, then normalize result
            var SuLength = Su.Length();
            if (!MathUtils.IsZero(SuLength))
            {
                Su /= SuLength;
            }

            var normal = Vector3.Cross(Su, Sv);
            var normalLength = normal.Length();

            return MathUtils.IsZero(normalLength) ? Vector3.Zero : normal / normalLength;
        }

        /// <inheritdoc cref="NurbsSurface{T}.GetCurvature(SurfaceCurvature, Vector2)"/>
        public static double ComputeCurvature(
            in NurbsSurface<Vector3> surface,
            SurfaceCurvature curvature,
            Vector2 uv
        )
        {
            var knotsU = surface.KnotsU;
            var knotsV = surface.KnotsV;

            Validate.Range(uv.X, knotsU[0], knotsU[knotsU.Count - 1], nameof(uv.X));
            Validate.Range(uv.Y, knotsV[0], knotsV[knotsV.Count - 1], nameof(uv.Y));

            // Get second-order rational derivatives
            var ders = ComputeRationalSurfaceDerivatives(in surface, 2, uv);

            var Su = ders[1][0]; // ∂S/∂u
            var Sv = ders[0][1]; // ∂S/∂v
            var Suu = ders[2][0]; // ∂²S/∂u²
            var Svv = ders[0][2]; // ∂²S/∂v²
            var Suv = ders[1][1]; // ∂²S/∂u∂v

            // Compute surface normal
            var normal = ComputeNormal(in surface, uv);

            // Second fundamental form coefficients
            double L = Vector3.Dot(Suu, normal);
            double M = Vector3.Dot(Suv, normal);
            double N = Vector3.Dot(Svv, normal);

            // First fundamental form coefficients
            double E = Vector3.Dot(Su, Su);
            double F = Vector3.Dot(Su, Sv);
            double G = Vector3.Dot(Sv, Sv);

            double denominator = E * G - F * F;
            if (MathUtils.IsZero(denominator))
            {
                return 0.0;
            }

            // Gaussian curvature
            double K = (L * N - M * M) / denominator;
            // Mean curvature
            double H = (E * N + G * L - 2 * F * M) / (2 * denominator);
            // Principal curvatures
            double discriminant = Math.Abs(H * H - K);
            double sqrtDisc = Math.Sqrt(discriminant);
            double k1 = H + sqrtDisc;
            double k2 = H - sqrtDisc;

            return curvature switch
            {
                SurfaceCurvature.Gauss => K,
                SurfaceCurvature.Mean => H,
                SurfaceCurvature.Maximum => k1,
                SurfaceCurvature.Minimum => k2,
                SurfaceCurvature.Abs => Math.Abs(k1) + Math.Abs(k2),
                SurfaceCurvature.Rms => Math.Sqrt(k1 * k1 + k2 * k2),
                _ => 0.0,
            };
        }

        /// <inheritdoc cref="NurbsSurface{T}.SwapDim()"/>
        public static NurbsSurface<Vector3> SwapDim(in NurbsSurface<Vector3> surface)
        {
            var controlPoints = surface.ControlPoints;
            int countU = controlPoints.Count;
            int countV = controlPoints[0].Count;

            // Transpose control points
            var transposedControlPoints = new ControlPoint<Vector3>[countV][];
            for (int j = 0; j < countV; j++)
            {
                transposedControlPoints[j] = new ControlPoint<Vector3>[countU];
                for (int i = 0; i < countU; i++)
                {
                    transposedControlPoints[j][i] = controlPoints[i][j];
                }
            }

            return new NurbsSurface<Vector3>(
                surface.DegreeV, // New DegreeU = old DegreeV
                surface.DegreeU, // New DegreeV = old DegreeU
                transposedControlPoints,
                surface.KnotsV, // New KnotsU = old KnotsV
                surface.KnotsU // New KnotsV = old KnotsU
            );
        }

        /// <inheritdoc cref="NurbsSurface{T}.Reverse(SurfaceDirection)"/>
        public static NurbsSurface<Vector3> Reverse(
            in NurbsSurface<Vector3> surface,
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
            var newControlPoints = new ControlPoint<Vector3>[countU][];

            if (direction == SurfaceDirection.UDirection)
            {
                // Reverse each column (swap rows for each column index)
                for (int i = 0; i < countU; i++)
                {
                    newControlPoints[i] = new ControlPoint<Vector3>[countV];
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
                    newControlPoints[i] = new ControlPoint<Vector3>[countV];
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
                    newControlPoints[i] = new ControlPoint<Vector3>[countV];
                    for (int j = 0; j < countV; j++)
                    {
                        newControlPoints[i][j] = controlPoints[countU - 1 - i][countV - 1 - j];
                    }
                }
            }

            return new NurbsSurface<Vector3>(
                surface.DegreeU,
                surface.DegreeV,
                newControlPoints,
                newKnotsU,
                newKnotsV
            );
        }

        /// <inheritdoc cref="NurbsSurface{T}.InsertKnot(double, int, SurfaceDirection, out NurbsSurface{T})"/>
        public static int InsertKnot(
            in NurbsSurface<Vector3> surface,
            double insertKnot,
            int times,
            SurfaceDirection direction,
            out NurbsSurface<Vector3> result
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

            var temp = new ControlPoint<Vector3>[degree + 1];
            int rows = controlPoints.Count;
            int columns = controlPoints[0].Count;

            if (isUDirection)
            {
                var updatedControlPoints = new ControlPoint<Vector3>[rows + times][];
                for (int i = 0; i < rows + times; i++)
                {
                    updatedControlPoints[i] = new ControlPoint<Vector3>[columns];
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

                result = new NurbsSurface<Vector3>(
                    surface.DegreeU,
                    surface.DegreeV,
                    updatedControlPoints,
                    insertedKnotVector,
                    surface.KnotsV
                );
            }
            else
            {
                var updatedControlPoints = new ControlPoint<Vector3>[rows][];
                for (int i = 0; i < rows; i++)
                {
                    updatedControlPoints[i] = new ControlPoint<Vector3>[columns + times];
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

                result = new NurbsSurface<Vector3>(
                    surface.DegreeU,
                    surface.DegreeV,
                    updatedControlPoints,
                    surface.KnotsU,
                    insertedKnotVector
                );
            }

            return times;
        }

        public static NurbsSurface<Vector3> RefineKnotVector(
            in NurbsSurface<Vector3> surface,
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
                var transposed = ControlPointsHelper.Transpose(controlPoints);
                var tempControlPoints = new List<ControlPoint<Vector3>[]>();
                IReadOnlyList<double> newKnotsU = null;

                for (int i = 0; i < transposed.Count; i++)
                {
                    var curve = new NurbsCurve<Vector3>(
                        surface.DegreeU,
                        transposed[i],
                        surface.KnotsU
                    );
                    var refined = NurbsCurve3DHelper.RefineKnotVector(curve, insertKnotElements);

                    var cpArray = new ControlPoint<Vector3>[refined.ControlPoints.Count];
                    for (int j = 0; j < refined.ControlPoints.Count; j++)
                    {
                        cpArray[j] = refined.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsU = refined.Knots;
                }

                // Convert to array and transpose back
                var tempArray = tempControlPoints.ToArray();
                var updatedControlPoints = ControlPointsHelper.Transpose(tempArray);

                return new NurbsSurface<Vector3>(
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
                var tempControlPoints = new List<ControlPoint<Vector3>[]>();
                IReadOnlyList<double> newKnotsV = null;

                for (int i = 0; i < controlPoints.Count; i++)
                {
                    var rowCPs = new ControlPoint<Vector3>[controlPoints[i].Count];
                    for (int j = 0; j < controlPoints[i].Count; j++)
                    {
                        rowCPs[j] = controlPoints[i][j];
                    }

                    var curve = new NurbsCurve<Vector3>(surface.DegreeV, rowCPs, surface.KnotsV);
                    var refined = NurbsCurve3DHelper.RefineKnotVector(curve, insertKnotElements);

                    var cpArray = new ControlPoint<Vector3>[refined.ControlPoints.Count];
                    for (int j = 0; j < refined.ControlPoints.Count; j++)
                    {
                        cpArray[j] = refined.ControlPoints[j];
                    }
                    tempControlPoints.Add(cpArray);
                    newKnotsV = refined.Knots;
                }

                return new NurbsSurface<Vector3>(
                    surface.DegreeU,
                    surface.DegreeV,
                    tempControlPoints.ToArray(),
                    surface.KnotsU,
                    newKnotsV
                );
            }
        }
    }
}
