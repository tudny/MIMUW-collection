package concurrentcube.sequantial;

import java.util.Arrays;

class Side {
    private final int[][] tiles;

    public Side(int size, int color) {
        this(size);

        for (int[] tile : tiles) {
            Arrays.fill(tile, color);
        }
    }

    public Side(int size) {
        this.tiles = new int[size][];

        for (int i = 0; i < size; ++i) {
            this.tiles[i] = new int[size];
        }
    }

    public String toStringInterrupted() throws InterruptedException {
        StringBuilder result = new StringBuilder();

        for (int y = this.tiles.length - 1; y >= 0 ; --y) {
            for (int x = 0; x < this.tiles.length; ++x) {
                if (Thread.currentThread().isInterrupted())
                    throw new InterruptedException();

                result.append(getTile(new Coordinates(x, y)));
            }
        }

        return result.toString();
    }

    public int getTile(Coordinates coordinates) {
        return tiles[coordinates.getX()][coordinates.getY()];
    }

    public void setTile(Coordinates coordinates, int color) {
        tiles[coordinates.getX()][coordinates.getY()] = color;
    }

    private int getSize() {
        return tiles.length;
    }

    public Side rotate() {
        return rotate(1);
    }

    public Side rotateCounter() {
        return rotateCounter(1);
    }

    private interface RotationInterface {
        Coordinates perform(Coordinates coordinates, int system, int howManyTimes);
    }

    public Side rotate(int howManyTimes) {
        return rotateAny(howManyTimes, Coordinates::performRotate);
    }

    public Side rotateCounter(int howManyTimes) {
        return rotateAny(howManyTimes, Coordinates::performRotateCounter);
    }

    private Side rotateAny(int howManyTimes, RotationInterface rotationInterface) {
        Side newSide = new Side(getSize());

        for (int x = 0; x < getSize(); ++x) {
            for (int y = 0; y < getSize(); ++y) {
                Coordinates current = new Coordinates(x, y);
                newSide.setTile(rotationInterface.perform(current, getSize(), howManyTimes), getTile(current));
            }
        }

        return newSide;
    }
}
