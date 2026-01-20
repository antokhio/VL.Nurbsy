using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    public static class NurbsSurface2DHelper
    {
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

        /// <summary>
        /// The NURBS Book 2nd Edition Page111
        /// Algorithm A3.6
        /// Compute surface derivatives.
        /// </summary>
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

        /// <summary>
        /// Optimized computation for first-order surface derivatives.
        /// </summary>
        /// <param name="surface">The NURBS surface.</param>
        /// <param name="uv">The uv parameter.</param>
        /// <returns>2x2 array of first-order derivatives.</returns>
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
    }
}
