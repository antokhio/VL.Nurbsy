using System.Collections.Immutable;
using System.Text.Json;
using Nurbsy;
using Nurbsy.Rendering.Helpers;
using Stride.Core.Mathematics;
using Xunit;

namespace Nurbsy.Tests;

public sealed class NurbsSurfaceSaveReloadTests
{
    private const float Epsilon = 1e-3f;
    private static readonly Int2 Tesselation = new(32, 32);
    private static readonly Vector2 Delta = new(0.3f, 1.5f);

    // Persisted form of a surface, mirrors what an editor would write to disk.
    private sealed record SavedSurface(
        int DegreeU,
        int DegreeV,
        float[][][] ControlPoints,
        double[] KnotsU,
        double[] KnotsV);

    // 3x3 quadratic, clamped knots, grid over (-8,-5)..(8,5). controlPoints[u][v].
    private static Vector2[][] Grid()
    {
        var grid = new Vector2[3][];
        for (int u = 0; u < 3; u++)
        {
            grid[u] = new Vector2[3];
            for (int v = 0; v < 3; v++)
                grid[u][v] = new Vector2(-8f + 8f * u, -5f + 5f * v);
        }
        return grid;
    }

    private static NurbsSurface<Vector2> Create(Vector2[][] grid, double[] knotsU, double[] knotsV)
    {
        var rows = grid
            .Select(col => (IReadOnlyList<Vector2>)col.ToImmutableArray())
            .ToImmutableArray();
        return new NurbsSurface<Vector2>(2, 2, rows, knotsU.ToImmutableArray(), knotsV.ToImmutableArray());
    }

    private static NurbsSurface<Vector2> Create(Vector2[][] grid)
    {
        double[] knots = [0, 0, 0, 1, 1, 1];
        return Create(grid, knots, knots);
    }

    // Horizontal middle row: v == 1 for every u.
    private static Vector2[][] MoveHorizontalMiddle(Vector2[][] grid, Vector2 delta)
    {
        var moved = grid.Select(col => (Vector2[])col.Clone()).ToArray();
        for (int u = 0; u < moved.Length; u++)
            moved[u][1] += delta;
        return moved;
    }

    private static string Save(NurbsSurface<Vector2> s) =>
        JsonSerializer.Serialize(new SavedSurface(
            s.DegreeU,
            s.DegreeV,
            s.ControlPoints
                .Select(col => col.Select(cp => new[] { cp.Value.X, cp.Value.Y, (float)cp.Weight }).ToArray())
                .ToArray(),
            s.KnotsU.ToArray(),
            s.KnotsV.ToArray()));

    private static NurbsSurface<Vector2> Load(string json)
    {
        var d = JsonSerializer.Deserialize<SavedSurface>(json)!;
        var rows = d.ControlPoints
            .Select(col => (IReadOnlyList<ControlPoint<Vector2>>)col
                .Select(p => new ControlPoint<Vector2>(new Vector2(p[0], p[1]), p[2]))
                .ToImmutableArray())
            .ToImmutableArray();
        return new NurbsSurface<Vector2>(
            d.DegreeU, d.DegreeV, rows, d.KnotsU.ToImmutableArray(), d.KnotsV.ToImmutableArray());
    }

    private static Vector3[] Mesh(NurbsSurface<Vector2> s, NurbsSurfaceArcLengthCache<Vector2> cache) =>
        NurbsSurfaceModelFactory.GenerateMeshDataCached(s, Tesselation, cache)
            .Vertices.Select(v => v.Position).ToArray();

    private static void AssertSameMesh(Vector3[] expected, Vector3[] actual, string context)
    {
        Assert.Equal(expected.Length, actual.Length);
        float worst = 0;
        int worstIndex = -1;
        for (int i = 0; i < expected.Length; i++)
        {
            var e = Vector3.Distance(expected[i], actual[i]);
            if (e > worst) { worst = e; worstIndex = i; }
        }
        Assert.True(worst < Epsilon, $"{context}: drift {worst} at vertex {worstIndex}");
    }

    [Fact]
    public void Save_and_load_round_trips_the_control_net()
    {
        var surface = Create(MoveHorizontalMiddle(Grid(), Delta));
        var reloaded = Load(Save(surface));

        for (int u = 0; u < 3; u++)
        for (int v = 0; v < 3; v++)
            Assert.Equal(surface.ControlPoints[u][v], reloaded.ControlPoints[u][v]);
        Assert.Equal(surface.KnotsU, reloaded.KnotsU);
        Assert.Equal(surface.KnotsV, reloaded.KnotsV);
    }

    [Fact]
    public void Moving_horizontal_middle_points_then_reload_keeps_the_shape()
    {
        var live = new NurbsSurfaceArcLengthCache<Vector2>();

        // Live session: rendered once untouched, then the middle row is moved.
        Mesh(Create(Grid()), live);
        var moved = Create(MoveHorizontalMiddle(Grid(), Delta));
        var liveMesh = Mesh(moved, live);

        // Save, then recreate with the same parameters in a fresh session.
        var reloaded = Load(Save(moved));
        var reloadedMesh = Mesh(reloaded, new NurbsSurfaceArcLengthCache<Vector2>());

        AssertSameMesh(liveMesh, reloadedMesh, "live vs reload");
    }

    [Fact]
    public void Incremental_drag_then_repeated_resave_keeps_the_shape()
    {
        var live = new NurbsSurfaceArcLengthCache<Vector2>();
        var grid = Grid();
        var surface = Create(grid);
        Mesh(surface, live);

        // Drag in small steps, rendering every frame with the same cache.
        const int steps = 10;
        for (int i = 0; i < steps; i++)
        {
            grid = MoveHorizontalMiddle(grid, Delta / steps);
            surface = Create(grid);
            Mesh(surface, live);
        }
        var liveMesh = Mesh(surface, live);

        // Save / reload several times; every reload must render the same shape.
        var json = Save(surface);
        for (int cycle = 0; cycle < 5; cycle++)
        {
            var reloaded = Load(json);
            AssertSameMesh(liveMesh, Mesh(reloaded, new NurbsSurfaceArcLengthCache<Vector2>()), $"resave {cycle}");
            json = Save(reloaded);
        }
    }

    [Fact]
    public void Unchanged_surface_reuses_cached_tables()
    {
        var cache = new NurbsSurfaceArcLengthCache<Vector2>();
        var surface = Create(Grid());

        cache.GetOrBuild(in surface);
        var u = cache.ULut;
        var v = cache.VLut;
        cache.GetOrBuild(in surface);

        Assert.Same(u, cache.ULut);
        Assert.Same(v, cache.VLut);
    }
}
