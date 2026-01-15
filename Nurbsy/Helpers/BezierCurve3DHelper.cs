using Nurbsy.Algorithm;
using Stride.Core.Mathematics;

namespace Nurbsy.Helpers
{
    internal static class BezierCurve3DHelper
    {
        public static Vector3 GetPointOnCurveByDeCasteljau(
            BezierCurve<Vector3> curve,
            double paramT
        )
        {
            int degree = curve.Degree;
            ControlPoint<Vector3>[] temp = curve.ControlPoints.ToArray();

            for (int k = 1; k <= degree; k++)
            {
                for (int i = 0; i <= degree - k; i++)
                {
                    var p0 = temp[i];
                    var p1 = temp[i + 1];

                    double w0 = p0.Weight;
                    double w1 = p1.Weight;

                    Vector3 v0 = p0.Value * (float)w0;
                    Vector3 v1 = p1.Value * (float)w1;

                    double t = paramT;
                    double oneMinusT = 1.0 - t;

                    double wNew = oneMinusT * w0 + t * w1;
                    Vector3 vNew = (float)oneMinusT * v0 + (float)t * v1;

                    if (Math.Abs(wNew) < MathUtils.Epsilon)
                    {
                        temp[i] = new ControlPoint<Vector3>(Vector3.Zero, (float)wNew);
                    }
                    else
                    {
                        temp[i] = new ControlPoint<Vector3>(vNew / (float)wNew, (float)wNew);
                    }
                }
            }

            return temp[0].Value;
        }

        public static Vector3 GetPointOnCurveByBernstein(BezierCurve<Vector3> curve, double paramT)
        {
            int degree = curve.Degree;
            var controlPoints = curve.ControlPoints;
            double[] bernsteinArray = Polynomials.AllBernstein(degree, paramT);

            double wDenominator = 0.0;
            Vector3 pNumerator = Vector3.Zero;

            for (int k = 0; k <= degree; k++)
            {
                var cp = controlPoints[k];
                double b = bernsteinArray[k];
                double w = cp.Weight;
                pNumerator += (float)(b * w) * cp.Value;
                wDenominator += b * w;
            }

            if (Math.Abs(wDenominator) > Constants.DoubleEpsilon)
            {
                pNumerator /= (float)wDenominator;
            }
            return pNumerator;
        }

        public static Vector3 GetPointOnQuadraticArc(
            ref ControlPoint<Vector3> startPoint,
            ref ControlPoint<Vector3> middlePoint,
            ref ControlPoint<Vector3> endPoint,
            double paramT
        )
        {
            Validate.Range(paramT, 0.0, 1.0, nameof(paramT));

            var p0 = startPoint.Value;
            var p1 = middlePoint.Value;
            var p2 = endPoint.Value;

            double w0 = startPoint.Weight;
            double w1 = middlePoint.Weight;
            double w2 = endPoint.Weight;

            double t = paramT;
            double t1 = 1.0 - t;
            double t1sq = t1 * t1;
            double tsq = t * t;
            double tt1 = 2 * t * t1;

            // Rational Bezier evaluation
            // P(t) = ( (1-t)^2 w0 P0 + 2t(1-t) w1 P1 + t^2 w2 P2 ) / ( (1-t)^2 w0 + 2t(1-t) w1 + t^2 w2 )

            double den = t1sq * w0 + tt1 * w1 + tsq * w2;
            if (Math.Abs(den) < MathUtils.Epsilon)
                return startPoint.Value;

            Vector3 num = (float)(t1sq * w0) * p0 + (float)(tt1 * w1) * p1 + (float)(tsq * w2) * p2;
            Vector3 res = num / (float)den;

            return res;
        }

        public static bool ComputeMiddleControlPointsOnQuadraticCurve(
            ref Vector3 startPoint,
            ref Vector3 startTangent,
            ref Vector3 endPoint,
            ref Vector3 endTangent,
            ref IList<ControlPoint<Vector3>> controlPoints
        )
        {
            Vector3 chord = (endPoint - startPoint);
            chord.Normalize();
            Vector3 nST = startTangent;
            nST.Normalize();
            Vector3 nET = endTangent;
            nET.Normalize();

            if (MathUtils.IsAlmostEqualTo(nST, chord) || MathUtils.IsAlmostEqualTo(nST, -chord))
            {
                if (MathUtils.IsAlmostEqualTo(nST, nET) || MathUtils.IsAlmostEqualTo(nST, -nET))
                {
                    controlPoints.Add(new ControlPoint<Vector3>((startPoint + endPoint) * 0.5f, 1));
                    return true;
                }
            }

            double alf1 = 0.0;
            double alf2 = 0.0;
            Vector3 R = Vector3.Zero;

            var type = Intersection3D.ComputeRays(
                startPoint,
                nST,
                endPoint,
                nET,
                out alf1,
                out alf2,
                out R
            );

            if (type == CurveCurveIntersectionType.Intersecting)
            {
                if (
                    MathUtils.IsAlmostEqualTo(startPoint, R)
                    || MathUtils.IsAlmostEqualTo(endPoint, R)
                )
                {
                    R = (startPoint + endPoint) * 0.5f;
                    double weight = 0.0;
                    if (
                        !ComputeWeightForRationalQuadraticInterpolation(
                            startPoint,
                            R,
                            endPoint,
                            out weight
                        )
                    )
                        return false;
                    controlPoints.Add(new ControlPoint<Vector3>(R, (float)weight));
                    return true;
                }

                if (MathUtils.IsGreaterThan(alf1, 0.0) && MathUtils.IsLessThan(alf2, 0.0))
                {
                    double weight = 0.0;
                    if (
                        !ComputeWeightForRationalQuadraticInterpolation(
                            startPoint,
                            R,
                            endPoint,
                            out weight
                        )
                    )
                        return false;
                    controlPoints.Add(new ControlPoint<Vector3>(R, (float)weight));
                    return true;
                }
            }

            Vector3 SE = endPoint - startPoint;
            Vector3 seNorm = SE;
            seNorm.Normalize();

            if (MathUtils.IsAlmostEqualTo(seNorm, nST) && MathUtils.IsAlmostEqualTo(seNorm, nET))
            {
                R = (startPoint + endPoint) * 0.5f;
                double weight = 0.0;
                if (
                    !ComputeWeightForRationalQuadraticInterpolation(
                        startPoint,
                        R,
                        endPoint,
                        out weight
                    )
                )
                    return false;
                controlPoints.Add(new ControlPoint<Vector3>(R, (float)weight));
                return true;
            }

            double gamma1 = 0.0;
            double gamma2 = 0.0;

            if (type != CurveCurveIntersectionType.Intersecting)
            {
                gamma1 = 0.5 * SE.Length();
                gamma2 = gamma1;
            }
            else
            {
                Vector3 SR = R - startPoint;
                Vector3 ER = R - endPoint;

                double theta0 = Vector3.Dot(SR, SE) / (double)(SR.Length() * SE.Length());
                double theta1 = Vector3.Dot(ER, SE) / (double)(ER.Length() * SE.Length());

                double alpha = 2.0 / 3.0;
                gamma1 = 0.5 * SE.Length() / (1.0 + alpha * theta1 + (1.0 - alpha) * theta0);
                gamma2 = 0.5 * SE.Length() / (1.0 + alpha * theta0 + (1.0 - alpha) * theta1);
            }

            Vector3 R1 = startPoint + (float)gamma1 * nST;
            Vector3 R2 = endPoint - (float)gamma2 * nET;
            Vector3 Qk = ((float)gamma1 * R2 + (float)gamma2 * R1) / (float)(gamma1 + gamma2);

            double weight1 = 0.0;
            if (!ComputeWeightForRationalQuadraticInterpolation(startPoint, R1, Qk, out weight1))
                return false;

            double weight2 = 0.0;
            if (!ComputeWeightForRationalQuadraticInterpolation(Qk, R2, endPoint, out weight2))
                return false;

            controlPoints.Add(new ControlPoint<Vector3>(R1, (float)weight1));
            controlPoints.Add(new ControlPoint<Vector3>(Qk, 1.0f));
            controlPoints.Add(new ControlPoint<Vector3>(R2, (float)weight2));

            return true;
        }

        private static bool ComputeWeightForRationalQuadraticInterpolation(
            Vector3 S,
            Vector3 R,
            Vector3 E,
            out double weight
        )
        {
            Vector3 chord = E - S;
            Vector3 tan = R - S;
            if (MathUtils.IsZero(chord) || MathUtils.IsZero(tan))
            {
                weight = 1.0;
                return true;
            }
            chord.Normalize();
            tan.Normalize();
            double dot = Vector3.Dot(chord, tan);

            double w = Math.Sqrt(Math.Max(0.0, 1.0 - dot * dot));

            weight = w > MathUtils.Epsilon ? w : 1.0;
            return true;
        }
    }
}
