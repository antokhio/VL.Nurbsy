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
        }

        public NurbsCurve(int degree, IReadOnlyList<T> controlPoints, IReadOnlyList<double> knots)
        {
            Degree = degree;
            ControlPoints = controlPoints.Select(cp => new ControlPoint<T>(cp)).ToArray();
            Knots = knots;
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
    }
}
