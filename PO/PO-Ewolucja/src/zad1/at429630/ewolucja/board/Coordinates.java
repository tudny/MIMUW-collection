package zad1.at429630.ewolucja.board;

import zad1.at429630.ewolucja.config.Config;

public class Coordinates {
    private final int x, y;
    private final Config config;

    public Coordinates(int x, int y, Config config) {
        this.x = x;
        this.y = y;
        this.config = config;
    }

    public Coordinates from(Direction... directions) {
        int newX = x;
        int newY = y;

        for (Direction direction : directions) {
            newX = (newX + direction.getDx() + config.getRozmiarPlanszyX()) % config.getRozmiarPlanszyX();
            newY = (newY + direction.getDy() + config.getRozmiarPlanszyY()) % config.getRozmiarPlanszyY();
        }
        return new Coordinates(newX, newY, config);
    }

    public int getX() {
        return x;
    }

    public int getY() {
        return y;
    }

    public static Coordinates randomCoordinates(Config config) {
        int x = config.getRandomizer().getRng(config.getRozmiarPlanszyX());
        int y = config.getRandomizer().getRng(config.getRozmiarPlanszyY());
        return new Coordinates(x, y, config);
    }

    @Override
    public String toString() {
        return "Coordinates{" +
                "x=" + x +
                ", y=" + y +
                '}';
    }
}
