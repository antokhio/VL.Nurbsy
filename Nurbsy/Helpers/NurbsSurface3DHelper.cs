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
    }
}
