namespace Nurbsy
{
    public enum CurveCurveIntersectionType
    {
        Intersecting = 0,
        Parallel = 1,
        Coincident = 2,
        Skew = 3,
    }

    public enum LinePlaneIntersectionType
    {
        Intersecting = 0,
        Parallel = 1,
        On = 2,
    }

    public enum CurveNormal
    {
        Normal = 0,
        Binormal = 1,
    }

    public enum SurfaceDirection
    {
        All = 0,
        UDirection = 1,
        VDirection = 2,
    }

    public enum SurfaceCurvature
    {
        Maximum = 0,
        Minimum = 1,
        Gauss = 2,
        Mean = 3,
        Abs = 4,
        Rms = 5,
    }

    public enum IntegratorType
    {
        Simpson = 0,
        GaussLegendre = 1,
        Chebyshev = 2,
    }

    public enum OffsetType
    {
        TillerAndHanson = 0,
        PieglAndTiller = 1,
    }
}
