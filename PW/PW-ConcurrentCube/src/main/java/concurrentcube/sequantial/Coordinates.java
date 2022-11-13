package concurrentcube.sequantial;

import static concurrentcube.sequantial.SequentialCube.ROTATIONS;

public class Coordinates {
    private final int x, y;

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }

    public Coordinates(int x, int y) {
        this.x = x;
        this.y = y;
    }

    private static int flip(int x, int system) {
        return system - x - 1;
    }

    private static int flipRotations(int howManyTimes) {
        return (((ROTATIONS - howManyTimes) % ROTATIONS) + ROTATIONS) % ROTATIONS;
    }

    public Coordinates rotate(int system) {
        return rotateSystemClockwise(this, system);
    }

    public Coordinates performRotate(int system, int howManyTime) {
        return performRotation(this, system, howManyTime);
    }

    public Coordinates performRotateCounter(int system, int howManyTime) {
        return performRotation(this, system, flipRotations(howManyTime));
    }

    public static Coordinates performRotation(Coordinates coordinates, int system, int howManyTimes) {
        while (howManyTimes --> 0)
            coordinates = rotateSystemClockwise(coordinates, system);

        return coordinates;
    }

    public static Coordinates rotateSystemClockwise(Coordinates coordinates, int system) {
        int newX = coordinates.getY();
        int newY = flip(coordinates.getX(), system);

        return new Coordinates(newX, newY);
    }

    public static Coordinates rotateSystemCounterClockwise(Coordinates coordinates, int system) {
        return performRotation(coordinates, system, flipRotations(1));
    }

    @Override
    public String toString() {
        return "(" +
                "" + x +
                ", " + y +
                ')';
    }
}
