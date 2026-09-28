using System.Collections.Immutable;
using System.Text.Json;
using Nurbsy.Algorithm;
using Nurbsy.Rendering.Helpers;
using Stride.Core.Mathematics;
using Xunit;

namespace Nurbsy.Tests;

public sealed class NurbsSurfaceMiddlePointDiagnostics
{
    private const float Epsilon = 1e-4f;
    private const int Tess = 33;
    private static readonly Int2 Tesselation = new(Tess, Tess);

    private sealed record SavedSurface(
        int DegreeU,
        int DegreeV,
        float[][][] ControlPoints,
        double[] KnotsU,
        double[] KnotsV
    );

    // n x n grid over (-0.5,-0.5)..(0.5,0.5), controlPoints[u][v].
    private static Vector3[][] Grid(int n) =>
        Enumerable
            .Range(0, n)
            .Select(u =>
                Enumerable
                    .Range(0, n)
                    .Select(v => new Vector3(u / (float)(n - 1) - 0.5f, v / (float)(n - 1) - 0.5f, 0))
                    .ToArray()
            )
            .ToArray();

    // Same construction path as NurbsSurfaceNode.Build(): immutable control points, clamped knots.
    private static NurbsSurface<Vector3> Create(int degree, Vector3[][] grid)
    {
        var rows = grid.Select(col =>
                (IReadOnlyList<ControlPoint<Vector3>>)
                    col.Select(p => new ControlPoint<Vector3>(p)).ToImmutableArray()
            )
            .ToImmutableArray();
        return new NurbsSurface<Vector3>(
            degree,
            degree,
            rows,
            KnotsUtils.GenerateClampedKnots(degree, grid.Length),
            KnotsUtils.GenerateClampedKnots(degree, grid[0].Length)
        );
    }

    private static Vector3[][] MoveMiddle(Vector3[][] grid, Vector3 delta)
    {
        var moved = grid.Select(col => (Vector3[])col.Clone()).ToArray();
        int m = grid.Length / 2;
        moved[m][m] += delta;
        return moved;
    }

    private static Vector3[] Mesh(NurbsSurface<Vector3> s) =>
        NurbsSurfaceModelFactory.GenerateMeshData(s, Tesselation).Vertices.Select(v => v.Position).ToArray();

    private static Vector3[] MeshCached(NurbsSurface<Vector3> s, NurbsSurfaceParamCache<Vector3> cache) =>
        NurbsSurfaceModelFactory.GenerateMeshDataCached(s, Tesselation, cache)
            .Vertices.Select(v => v.Position).ToArray();

    // Vertex at the texture coordinate (0.5, 0.5).
    private static Vector3 TextureCenter(Vector3[] mesh) => mesh[(Tess / 2) * Tess + Tess / 2];

    private static void AssertSameMesh(Vector3[] expected, Vector3[] actual, string context)
    {
        Assert.Equal(expected.Length, actual.Length);
        float worst = 0;
        int worstIndex = -1;
        for (int i = 0; i < expected.Length; i++)
        {
            var e = Vector3.Distance(expected[i], actual[i]);
            if (e > worst)
            {
                worst = e;
                worstIndex = i;
            }
        }
        Assert.True(worst < Epsilon, $"{context}: drift {worst} at vertex {worstIndex}");
    }

    private static string Save(NurbsSurface<Vector3> s) =>
        JsonSerializer.Serialize(
            new SavedSurface(
                s.DegreeU,
                s.DegreeV,
                s.ControlPoints.Select(col =>
                        col.Select(cp => new[] { cp.Value.X, cp.Value.Y, cp.Value.Z, (float)cp.Weight }).ToArray()
                    )
                    .ToArray(),
                s.KnotsU.ToArray(),
                s.KnotsV.ToArray()
            )
        );

    private static NurbsSurface<Vector3> Load(string json)
    {
        var d = JsonSerializer.Deserialize<SavedSurface>(json)!;
        var rows = d.ControlPoints.Select(col =>
                (IReadOnlyList<ControlPoint<Vector3>>)
                    col.Select(p => new ControlPoint<Vector3>(new Vector3(p[0], p[1], p[2]), p[3]))
                        .ToImmutableArray()
            )
            .ToImmutableArray();
        return new NurbsSurface<Vector3>(
            d.DegreeU,
            d.DegreeV,
            rows,
            d.KnotsU.ToImmutableArray(),
            d.KnotsV.ToImmutableArray()
        );
    }

    // Texture (0.5,0.5) must shift by N_mid(0.5)^2 * delta: the middle control point's basis weight
    // there is 1 for degree 1 and 0.5 * 0.5 for quadratic. Arc-length reparametrization cancels it.
    [Theory]
    [InlineData(1, 1f)]
    [InlineData(2, 0.25f)]
    public void Moving_middle_point_moves_the_texture_center(int degree, float weight)
    {
        var delta = new Vector3(0.3f, 0.2f, 0);
        var before = TextureCenter(Mesh(Create(degree, Grid(3))));
        var after = TextureCenter(Mesh(Create(degree, MoveMiddle(Grid(3), delta))));

        var shift = after - before;
        Assert.True(
            Vector3.Distance(shift, delta * weight) < Epsilon,
            $"texture center shift {shift}, expected {delta * weight}"
        );
    }

    [Theory]
    [InlineData(2, 3)]
    [InlineData(3, 4)]
    [InlineData(3, 6)]
    public void Uniform_grid_has_no_texture_compression(int degree, int count)
    {
        var mesh = Mesh(Create(degree, Grid(count)));
        for (int i = 0; i < Tess; i++)
        for (int j = 0; j < Tess; j++)
        {
            var expected = new Vector3(i / (float)(Tess - 1) - 0.5f, j / (float)(Tess - 1) - 0.5f, 0);
            Assert.True(
                Vector3.Distance(expected, mesh[j * Tess + i]) < Epsilon,
                $"vertex ({i},{j}) = {mesh[j * Tess + i]}, expected {expected}"
            );
        }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void Cached_mesh_matches_uncached_after_moving_middle(int degree)
    {
        int n = degree + 2;
        var cache = new NurbsSurfaceParamCache<Vector3>();
        MeshCached(Create(degree, Grid(n)), cache);

        var moved = Create(degree, MoveMiddle(Grid(n), new Vector3(0.3f, 0.2f, 0.1f)));
        AssertSameMesh(Mesh(moved), MeshCached(moved, cache), "cached vs uncached");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void Rebuilding_from_saved_params_reproduces_the_live_mesh(int degree)
    {
        int n = degree + 2;
        var live = new NurbsSurfaceParamCache<Vector3>();

        // Drag the middle point in small steps, rendering every frame with the same cache.
        var grid = Grid(n);
        var surface = Create(degree, grid);
        MeshCached(surface, live);
        for (int i = 0; i < 10; i++)
        {
            grid = MoveMiddle(grid, new Vector3(0.03f, 0.02f, 0.01f));
            surface = Create(degree, grid);
            MeshCached(surface, live);
        }
        var liveMesh = MeshCached(surface, live);

        // Rebuild from control points + degrees + knots only, with no cache.
        var json = Save(surface);
        for (int cycle = 0; cycle < 3; cycle++)
        {
            var reloaded = Load(json);
            for (int u = 0; u < n; u++)
            for (int v = 0; v < n; v++)
                Assert.Equal(surface.ControlPoints[u][v], reloaded.ControlPoints[u][v]);
            Assert.Equal(surface.KnotsU, reloaded.KnotsU);
            Assert.Equal(surface.KnotsV, reloaded.KnotsV);

            AssertSameMesh(liveMesh, Mesh(reloaded), $"reload {cycle}");
            json = Save(reloaded);
        }
    }

    [Fact]
    public void Param_cache_ignores_control_point_positions()
    {
        var cache = new NurbsSurfaceParamCache<Vector3>();
        cache.GetOrBuild(Create(2, Grid(3)), Tess, Tess);
        var u = cache.UParams;
        var v = cache.VParams;

        // New surface with new knot instances but equal topology still rebuilds (reference equality),
        // yet produces identical parameters regardless of the moved middle point.
        cache.GetOrBuild(Create(2, MoveMiddle(Grid(3), new Vector3(0.3f, 0.2f, 0))), Tess, Tess);
        Assert.Equal(u, cache.UParams);
        Assert.Equal(v, cache.VParams);
    }
}