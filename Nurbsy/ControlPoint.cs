namespace Nurbsy
{
    public record struct ControlPoint<T>
        where T : struct
    {
        public T Value { get; set; }
        public double Weight { get; set; }

        public ControlPoint(T value)
        {
            Value = value;
            Weight = 1f;
        }

        public ControlPoint(T value, float weight)
        {
            Value = value;
            Weight = weight;
        }
    }
}
