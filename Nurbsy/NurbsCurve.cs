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

        public NurbsCurve(
            int degree,
            IReadOnlyList<ControlPoint<T>> controlPoints,
            IReadOnlyList<double> knots
        )
        {
            Degree = degree;
            ControlPoints = controlPoints;
            Knots = knots;

            Check();
        }

        public NurbsCurve(int degree, IReadOnlyList<T> controlPoints, IReadOnlyList<double> knots)
        {
            Degree = degree;
            ControlPoints = controlPoints.Select(cp => new ControlPoint<T>(cp)).ToImmutableArray();
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

        public T GetPointOnCurveByCornerCut(double paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.GetPointOnCurveByCornerCut(curve2, paramT);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.GetPointOnCurveByCornerCut(curve3, paramT);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public IReadOnlyList<T> GetDerivatives(int derivative, double paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.ComputeRationalCurveDerivatives(
                    curve2,
                    derivative,
                    paramT
                );
                return Unsafe.As<IReadOnlyList<Vector2>, IReadOnlyList<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.ComputeRationalCurveDerivatives(
                    curve3,
                    derivative,
                    paramT
                );
                return Unsafe.As<IReadOnlyList<Vector3>, IReadOnlyList<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public bool CanComputeDerivative(double paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.CanComputeDerivative(curve2, paramT);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.CanComputeDerivative(curve3, paramT);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public double GetCurvature(double paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.GetCurvature(curve2, paramT);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.GetCurvature(curve3, paramT);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public double GetTorsion(double paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.GetTorsion(curve2, paramT);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.GetTorsion(curve3, paramT);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public T GetNormal(CurveNormal normalType, double paramT)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.GetNormal(curve2, normalType, paramT);
                return Unsafe.As<Vector2, T>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.GetNormal(curve3, normalType, paramT);
                return Unsafe.As<Vector3, T>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public IReadOnlyList<T> ProjectNormal()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.ProjectNormal(curve2);
                return Unsafe.As<IReadOnlyList<Vector2>, IReadOnlyList<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.ProjectNormal(curve3);
                return Unsafe.As<IReadOnlyList<Vector3>, IReadOnlyList<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public IReadOnlyList<BezierCurve<T>> DecomposeToBeziers()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.DecomposeToBeziers(curve2);
                return Unsafe.As<
                    IReadOnlyList<BezierCurve<Vector2>>,
                    IReadOnlyList<BezierCurve<T>>
                >(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.DecomposeToBeziers(curve3);
                return Unsafe.As<
                    IReadOnlyList<BezierCurve<Vector3>>,
                    IReadOnlyList<BezierCurve<T>>
                >(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public (
            IReadOnlyList<T> TessellatedPoints,
            IReadOnlyList<double> CorrespondingKnots
        ) EquallyTessellate()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var (pts, knots) = NurbsCurve2DHelper.EquallyTessellate(curve2);
                var ptsT = Unsafe.As<IReadOnlyList<Vector2>, IReadOnlyList<T>>(ref pts);
                return (ptsT, knots);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var (pts, knots) = NurbsCurve3DHelper.EquallyTessellate(curve3);
                var ptsT = Unsafe.As<IReadOnlyList<Vector3>, IReadOnlyList<T>>(ref pts);
                return (ptsT, knots);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public double GetParamOnCurve(T point)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                ref var p = ref Unsafe.As<T, Vector2>(ref point);
                return NurbsCurve2DHelper.GetParamOnCurve(curve2, p);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                ref var p = ref Unsafe.As<T, Vector3>(ref point);
                return NurbsCurve3DHelper.GetParamOnCurve(curve3, p);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public NurbsCurve<T> Reparametrize(double alpha, double beta, double gamma, double delta)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.Reparametrize(curve2, alpha, beta, gamma, delta);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.Reparametrize(curve3, alpha, beta, gamma, delta);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public NurbsCurve<T> Reparametrize(double min, double max)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.Reparametrize(curve2, min, max);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.Reparametrize(curve3, min, max);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public NurbsCurve<T> RefineKnotVector(IReadOnlyList<double> insertKnotElements)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var refined2 = NurbsCurve2DHelper.RefineKnotVector(curve2, insertKnotElements);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref refined2);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var refined3 = NurbsCurve3DHelper.RefineKnotVector(curve3, insertKnotElements);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref refined3);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public NurbsCurve<T> Reverse()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.Reverse(curve2);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.Reverse(curve3);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public bool SplitAt(double parameter, out NurbsCurve<T> left, out NurbsCurve<T> right)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.SplitAt(curve2, parameter, out var l2, out var r2))
                {
                    left = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref l2);
                    right = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref r2);
                    return true;
                }
                left = default;
                right = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.SplitAt(curve3, parameter, out var l3, out var r3))
                {
                    left = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref l3);
                    right = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref r3);
                    return true;
                }
                left = default;
                right = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public bool Segment(double startParameter, double endParameter, out NurbsCurve<T> segment)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.Segment(curve2, startParameter, endParameter, out var res2))
                {
                    segment = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                segment = default;
                return false;
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.Segment(curve3, startParameter, endParameter, out var res3))
                {
                    segment = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                segment = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public bool RemoveKnot(double removeKnot, int times, out NurbsCurve<T> result)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.RemoveKnot(curve2, removeKnot, times, out var res2))
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default; // Should be curve copy?
                return false;
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.RemoveKnot(curve3, removeKnot, times, out var res3))
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public NurbsCurve<T> RemoveExcessiveKnots()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.RemoveExcessiveKnots(curve2);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.RemoveExcessiveKnots(curve3);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public NurbsCurve<T> ElevateDegree(int times)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                var result = NurbsCurve2DHelper.ElevateDegree(curve2, times);
                return Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref result);
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                var result = NurbsCurve3DHelper.ElevateDegree(curve3, times);
                return Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref result);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public bool ReduceDegree(out NurbsCurve<T> result)
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                if (NurbsCurve2DHelper.ReduceDegree(curve2, out var res2))
                {
                    result = Unsafe.As<NurbsCurve<Vector2>, NurbsCurve<T>>(ref res2);
                    return true;
                }
                result = default;
                return false;
            }
            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                if (NurbsCurve3DHelper.ReduceDegree(curve3, out var res3))
                {
                    result = Unsafe.As<NurbsCurve<Vector3>, NurbsCurve<T>>(ref res3);
                    return true;
                }
                result = default;
                return false;
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }

        public bool IsLinear()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.IsLinear(curve2);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.IsLinear(curve3);
            }

            return false;
        }

        public bool IsClosed()
        {
            if (typeof(T) == typeof(Vector2))
            {
                ref var curve2 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector2>>(ref this);
                return NurbsCurve2DHelper.IsClosed(curve2);
            }

            if (typeof(T) == typeof(Vector3))
            {
                ref var curve3 = ref Unsafe.As<NurbsCurve<T>, NurbsCurve<Vector3>>(ref this);
                return NurbsCurve3DHelper.IsClosed(curve3);
            }

            throw new NotSupportedException($"Type {typeof(T).Name} is not supported.");
        }
    }
}
