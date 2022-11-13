package concurrentcube.sequantial;

public class SequentialCube {
    protected static final int SIDES = 6;
    public static final int ROTATIONS = 4;
    public static final SideInfo[] RELATION = {
            new SideInfo(5, new int[]{ 4, 3, 2, 1 }, new int[]{ 2, 2, 2, 2 }), // 0
            new SideInfo(3, new int[]{ 0, 2, 5, 4 }, new int[]{ 1, 1, 1, 3 }), // 1
            new SideInfo(4, new int[]{ 0, 3, 5, 1 }, new int[]{ 0, 1, 2, 3 }), // 2
            new SideInfo(1, new int[]{ 0, 4, 5, 2 }, new int[]{ 3, 1, 3, 3 }), // 3
            new SideInfo(2, new int[]{ 0, 1, 5, 3 }, new int[]{ 2, 1, 0, 3 }), // 4
            new SideInfo(0, new int[]{ 1, 2, 3, 4 }, new int[]{ 0, 0, 0, 0 }), // 5
    };

    private final Side[] tile = new Side[SIDES];
    private final int size;

    private int getSize() {
        return size;
    }

    public SequentialCube(int size) {
        this.size = size;
        for (int i = 0; i < SIDES; ++i) {
            this.tile[i] = new Side(size, i);
        }
    }

    public void rotate(int side, int layer) {
        // We read information about this face of the cube.
        final SideInfo sideInfo = RELATION[side];

        if (layer == 0) { // If the layer is also connected to the face we rotate it.
            tile[side] = tile[side].rotate();
        }

        if (layer + 1 == getSize()) { // Symmetrically the opposite face.
            tile[sideInfo.getOpposite()] = tile[sideInfo.getOpposite()].rotateCounter();
        }

        for (int x = 0; x < getSize(); ++x) { // For every tile of the layer we do a cyclic swap of every four tails.
            Coordinates current = new Coordinates(x, layer);
            Coordinates corresponding = current.performRotate(getSize(), sideInfo.getHowManyRotations()[ROTATIONS - 1]);
            int toPut = tile[sideInfo.getNear()[ROTATIONS - 1]]
                    .getTile(corresponding);

            for (int k = 0; k < ROTATIONS; ++k) { // Cyclic swap.
                corresponding = current.performRotate(getSize(), sideInfo.getHowManyRotations()[k]);
                int newToPut = tile[sideInfo.getNear()[k]].getTile(corresponding);
                tile[sideInfo.getNear()[k]].setTile(corresponding, toPut);
                toPut = newToPut;
            }
        }
    }

    public String show() throws InterruptedException {
        StringBuilder result = new StringBuilder();
        for (int side = 0; side < SIDES; ++side) {
            result.append(tile[side].toStringInterrupted());
        }
        return result.toString();
    }

    public static void main(String[] args) throws InterruptedException {
        SequentialCube cube = new SequentialCube(4);

        cube.rotate(2, 0);
        cube.rotate(5, 1);

        System.out.println(cube.show());
    }
}
