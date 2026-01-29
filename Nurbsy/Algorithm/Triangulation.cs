namespace Nurbsy.Algorithm
{
    public static class Triangulation
    {
        public static int[] GenerateGridIndeciesCW(int tesselationX, int tesselationY)
        {
            int count = (tesselationX - 1) * (tesselationY - 1) * 6;
            int[] indices = GC.AllocateUninitializedArray<int>(count);

            int index = 0;
            for (int y = 0; y < tesselationY - 1; y++)
            {
                int rowOffset = y * tesselationX;
                int nextRowOffset = (y + 1) * tesselationX;

                for (int x = 0; x < tesselationX - 1; x++)
                {
                    var bl = rowOffset + x; // Bottom-Left
                    var br = rowOffset + x + 1; // Bottom-Right
                    var tl = nextRowOffset + x; // Top-Left
                    var tr = nextRowOffset + x + 1; // Top-Right

                    // Triangle 1: Bottom-Left -> Top-Left -> Bottom-Right (CW)
                    indices[index++] = bl;
                    indices[index++] = tl;
                    indices[index++] = br;

                    // Triangle 2: Bottom-Right -> Top-Left -> Top-Right (CW)
                    indices[index++] = br;
                    indices[index++] = tl;
                    indices[index++] = tr;
                }
            }
            return indices;
        }

        public static int[] GenerateTapeIndeciesCW(int tesselation)
        {
            int count = tesselation * 6;
            int[] indices = GC.AllocateUninitializedArray<int>(count);
            int index = 0;
            for (int i = 0; i < tesselation; i++)
            {
                var bl = i * 2; // Bottom-Left
                var br = i * 2 + 1; // Bottom-Right
                var tl = (i + 1) * 2; // Top-Left
                var tr = (i + 1) * 2 + 1; // Top-Right
                // Triangle 1: Bottom-Left -> Top-Left -> Bottom-Right (CW)
                indices[index++] = bl;
                indices[index++] = tl;
                indices[index++] = br;
                // Triangle 2: Bottom-Right -> Top-Left -> Top-Right (CW)
                indices[index++] = br;
                indices[index++] = tl;
                indices[index++] = tr;
            }
            return indices;
        }
    }
}
