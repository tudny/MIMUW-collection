package zad1.at429630.ewolucja.board;

import zad1.at429630.ewolucja.utils.Randomizer;

public enum Direction {
    NORTH(0, 0, 1),
    EAST(1, 1, 0),
    SOUTH(2, 0, -1),
    WEST(3, -1, 0);

    private final int direction;
    private final int dx;
    private final int dy;

    Direction(int direction, int dx, int dy) {
        this.direction = direction;
        this.dx = dx;
        this.dy = dy;
    }

    public Direction next() {
        return of((direction + 1) % getDirectionCount());
    }

    public Direction prev() {
        return of((direction + getDirectionCount() - 1) % getDirectionCount());
    }

    public int getDx() {
        return dx;
    }

    public int getDy() {
        return dy;
    }

    public static int getDirectionCount() {
        return Direction.values().length;
    }

    public static Direction of(int x) {
        for (Direction direction : Direction.values()) {
            if (direction.direction == x) {
                return direction;
            }
        }

        throw new IllegalArgumentException("No direction of " + x + ".");
    }

    public static Direction of(int dx, int dy) {
        for (Direction direction : Direction.values()) {
            if (direction.dx == dx && direction.dy == dy) {
                return direction;
            }
        }

        return null;
    }

    public static Direction random(Randomizer randomizer) {
        int index = randomizer.getRng(4);
        return values()[index];
    }
}
