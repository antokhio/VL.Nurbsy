using System.Collections.Immutable;
using System.Runtime.CompilerServices;
using Nurbsy.Helpers;
using Stride.Core.Mathematics;

namespace Nurbsy
{
    public record struct NurbsCurve<T>
        where T : struct
    {
        public int Degree { get; }
        public IReadOnlyList<ControlPoint<T>> ControlPoints { get; }
        public IReadOnlyList<double> Knots { get; }

        public NurbsCurve(int degree, IReadOnlyList<ControlPoint<T>> controlPoints)
        {
            Degree = degree;
            ControlPoints = controlPoints;
            Knots = NurbsCurveHelper.GenerateClampedKnots(Degree, ControlPoints.Count);

            Check();
        }

        public NurbsCurve(int degree, IReadOnlyList<T> points)
        {
            Degree = degree;
            ControlPoints = points.Select(cp => new ControlPoint<T>(cp)).ToImmutableArray();
            Knots = NurbsCurveHelper.GenerateClampedKnots(Degree, ControlPoints.Count);

            Check();
        }

        public NurbsCurve(int degree, IReadOnlyList<T> points, IReadOnlyList<double> knots)
        {
            Degree = degree;
            ControlPoints = points.Select(cp => new ControlPoint<T>(cp)).ToImmutableArray();
            Knots = knots;

            Check();
        }

        public void Check()
        {
            Validate.Argument(Degree > 0, nameof(Degree), "Degree must be gretater then zero.");
            Validate.Argument(
                Knots.Count > 0,
                nameof(Knots),
                "Knot vector must contain at least one knot."
            );
            Validate.Argument(
                Validate.IsValidKnots(Knots),
                nameof(Knots),
                "Knot must be a nondecreasing sequence of real numbers."
            );
            Validate.Argument(
                ControlPoints.Count > 0,
                nameof(ControlPoints),
                "ControlPoints must contain one point at least."
            );
            Validate.Argument(
                Validate.IsValidNURBS(Degree, ControlPoints.Count, Knots.Count),
                nameof(NurbsCurve<T>),
                "Arguments must be fit: m = n + p + 1"
            );
        }

        public T GetPointOnCurve(double paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.GetPointOnCurve(curve2, paramT);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.GetPointOnCurve(curve3, paramT);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public bool IsLinear()
        {
            var count = ControlPoints.Count;
            if (count < 2)
                return false;
            if (count == 2)
                return true;

            if (typeof(T) == typeof(Vector2))
            {
                var p0_t = ControlPoints[0].Value;
                var p0 = Unsafe.As<T, Vector2>(ref p0_t);
                var v0 = Vector2.Zero;

                int i = 1;
                for (; i < ControlPoints.Count; i++)
                {
                    var pi_t = ControlPoints[i].Value;
                    var pi = Unsafe.As<T, Vector2>(ref pi_t);
                    v0 = pi - p0;
                    if (v0.LengthSquared() > MathUtil.ZeroTolerance)
                        break;
                }

                if (i == ControlPoints.Count)
                    return true;

                v0 = Vector2.Normalize(v0);

                for (int k = 2; k < ControlPoints.Count; k++)
                {
                    var pk_t = ControlPoints[k].Value;
                    var pk = Unsafe.As<T, Vector2>(ref pk_t);
                    var vk = pk - p0;

                    // 2D Cross Product (Z component)
                    float cross = v0.X * vk.Y - v0.Y * vk.X;
                    if (Math.Abs(cross) > MathUtil.ZeroTolerance)
                        return false;
                }
                return true;
            }
            if (typeof(T) == typeof(Vector3))
            {
                var p0_t = ControlPoints[0].Value;
                var p0 = Unsafe.As<T, Vector3>(ref p0_t);
                var v0 = Vector3.Zero;

                int i = 1;
                for (; i < ControlPoints.Count; i++)
                {
                    var pi_t = ControlPoints[i].Value;
                    var pi = Unsafe.As<T, Vector3>(ref pi_t);
                    v0 = pi - p0;
                    if (v0.LengthSquared() > MathUtil.ZeroTolerance)
                        break;
                }

                if (i == ControlPoints.Count)
                    return true;

                v0 = Vector3.Normalize(v0);

                for (int k = 2; k < ControlPoints.Count; k++)
                {
                    var pk_t = ControlPoints[k].Value;
                    var pk = Unsafe.As<T, Vector3>(ref pk_t);
                    var vk = pk - p0;

                    var cross = Vector3.Cross(v0, vk);
                    if (cross.LengthSquared() > MathUtil.ZeroTolerance * MathUtil.ZeroTolerance)
                        return false;
                }
                return true;
            }

            return false;
        }
    }
}
