using Nurbsy.Algorithm;
using Stride.Core.Mathematics;
using Stride.Graphics;

namespace Nurbsy.Rendering.Helpers
{
    public static class NurbsCurveTapeModelFactory
    {
        public static GeometricMeshData<VertexPositionNormalTexture> GenerateTape(
            NurbsCurve<Vector2> nurbsCurve,
            float width,
            int tesselation
        )
        {
            var tes = Math.Max(2, tesselation);

            var vertexCount = (tes + 1) * 2;
            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(vertexCount);

            var knots = nurbsCurve.Knots;
            double tStart = knots[nurbsCurve.Degree];
            double tEnd = knots[knots.Count - 1 - nurbsCurve.Degree];
            double tRange = tEnd - tStart;

            double halfWidth = width * 0.5;
            var normal = Vector3.UnitZ;

            // Generate vertices
            for (int i = 0; i <= tes; i++)
            {
                double ratio = (double)i / tes;
                double t = tStart + tRange * ratio;

                // Evaluate curve for Position and Tangent (Index 0: Point, Index 1: First Derivative)
                var derivatives = nurbsCurve.GetDerivatives(1, (float)t);
                var pos2 = derivatives[0];
                var tan2 = derivatives[1];

                // Handle zero-length tangent
                if (tan2.LengthSquared() < 1e-6f)
                {
                    if (i > 0)
                    {
                        var prevPos = vertices[(i - 1) * 2].Position;
                        tan2 = pos2 - new Vector2(prevPos.X, prevPos.Y);
                    }

                    if (tan2.LengthSquared() < 1e-6f)
                        tan2 = Vector2.UnitX;
                }

                tan2 = Vector2.Normalize(tan2);

                // Calculate binormal (width direction).
                // In 3D we do Cross(tangent, dir). Here dir is UnitZ.
                // Cross((tx, ty, 0), (0, 0, 1)) -> (ty, -tx, 0)
                var binormal = new Vector2(tan2.Y, -tan2.X);

                var pLeft2 = pos2 - binormal * (float)halfWidth;
                var pRight2 = pos2 + binormal * (float)halfWidth;

                var pLeft = new Vector3(pLeft2.X, pLeft2.Y, 0f);
                var pRight = new Vector3(pRight2.X, pRight2.Y, 0f);

                // Texture coordinates: U varies along width (0..1), V varies along length (0..1)
                vertices[i * 2] = new VertexPositionNormalTexture(
                    pLeft,
                    normal,
                    new Vector2(0, (float)ratio)
                );
                vertices[i * 2 + 1] = new VertexPositionNormalTexture(
                    pRight,
                    normal,
                    new Vector2(1, (float)ratio)
                );
            }

            // Generate Indices (CW winding)
            var indices = Triangulation.GenerateTapeIndeciesCW(tes);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }

        public static GeometricMeshData<VertexPositionNormalTexture> GenerateTape(
            NurbsCurve<Vector3> nurbsCurve,
            float width,
            int tesselation,
            Vector3? direction = null
        )
        {
            var dir = direction ?? Vector3.UnitZ;

            var vertexCount = (tesselation + 1) * 2;
            var vertices = GC.AllocateUninitializedArray<VertexPositionNormalTexture>(vertexCount);

            var knots = nurbsCurve.Knots;
            double tStart = knots[nurbsCurve.Degree];
            double tEnd = knots[knots.Count - 1 - nurbsCurve.Degree];
            double tRange = tEnd - tStart;

            double halfWidth = width * 0.5;

            // Store previous binormal for parallel transport
            Vector3 prevBinormal = Vector3.Zero;

            // Generate vertices
            for (int i = 0; i <= tesselation; i++)
            {
                double ratio = (double)i / tesselation;
                double t = tStart + tRange * ratio;

                // Evaluate curve for Position and Tangent (Index 0: Point, Index 1: First Derivative)
                var derivatives = nurbsCurve.GetDerivatives(1, (float)t);
                var position = derivatives[0];
                var tangent = derivatives[1];

                // Handle zero-length tangent
                if (tangent.LengthSquared() < 1e-6f)
                {
                    if (i > 0)
                        tangent = position - vertices[(i - 1) * 2].Position;

                    if (tangent.LengthSquared() < 1e-6f)
                        tangent = Vector3.UnitX;
                }

                tangent = Vector3.Normalize(tangent);

                Vector3 binormal;

                if (i == 0)
                {
                    // Initial frame calculation using reference direction
                    binormal = Vector3.Cross(tangent, dir);

                    // Handle parallel case
                    if (binormal.LengthSquared() < 1e-6f)
                    {
                        if (Math.Abs(tangent.Y) < 0.99f)
                            binormal = Vector3.Cross(tangent, Vector3.UnitY);
                        else
                            binormal = Vector3.Cross(tangent, Vector3.UnitX);
                    }
                    binormal = Vector3.Normalize(binormal);
                }
                else
                {
                    // Parallel Transport (Rotation Minimizing Frame)
                    // Instead of recalculating Cross(tangent, dir) which flips when tangent crosses dir,
                    // we project the previous binormal onto the new normal plane.
                    // This is similar to the logic in NurbsCurve3DHelper.ProjectNormal.

                    var dot = Vector3.Dot(prevBinormal, tangent);
                    binormal = prevBinormal - dot * tangent;

                    // If projection is too small (e.g. 90 degree sharp turn), fallback to reference
                    if (binormal.LengthSquared() < 1e-6f)
                    {
                        binormal = Vector3.Cross(tangent, dir);
                        if (binormal.LengthSquared() < 1e-6f)
                        {
                            if (Math.Abs(tangent.Y) < 0.99f)
                                binormal = Vector3.Cross(tangent, Vector3.UnitY);
                            else
                                binormal = Vector3.Cross(tangent, Vector3.UnitX);
                        }
                    }

                    binormal = Vector3.Normalize(binormal);
                }

                prevBinormal = binormal;
                var normal = Vector3.Cross(binormal, tangent);

                var pLeft = position - binormal * (float)halfWidth;
                var pRight = position + binormal * (float)halfWidth;

                // Texture coordinates: U varies along width (0..1), V varies along length (0..1)
                vertices[i * 2] = new VertexPositionNormalTexture(
                    pLeft,
                    normal,
                    new Vector2(0, (float)ratio)
                );
                vertices[i * 2 + 1] = new VertexPositionNormalTexture(
                    pRight,
                    normal,
                    new Vector2(1, (float)ratio)
                );
            }

            // Generate Indices (CW winding)
            var indices = Triangulation.GenerateTapeIndeciesCW(tesselation);

            return new GeometricMeshData<VertexPositionNormalTexture>(
                vertices,
                indices,
                isLeftHanded: false
            );
        }
    }
}
