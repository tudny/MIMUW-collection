package zad1.at429630.ewolucja.board;

import zad1.at429630.ewolucja.board.field.Field;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class Board {
    private final Field[][] fields;
    private final int x;
    private final int y;

    public Board(int x, int y, List<String> fieldFileAsLines) {
        if (x <= 0 || y <= 0)
            throw new IllegalArgumentException(
                    "Board shouldn't be degenerated. The measurements are x:" + x + " y:" + y + ".");

        this.x = x;
        this.y = y;
        fields = new Field[x][y];

        int actY = 0;

        for (String line : fieldFileAsLines) {
            if (line.length() != x)
                throw new IllegalArgumentException("Line lengths of lines should be the same.");

            int actX = 0;

            for (Character character : line.toCharArray()) {
                setField(character, actX, actY);
                actX++;
            }
            actY++;
        }
    }

    private void setField(Character character, int x, int y) {
        fields[x][y] = Field.of(character);
    }

    public Field getField(Coordinates coordinates) {
        return fields[coordinates.getX()][coordinates.getY()];
    }

    public int getWidth() {
        return x;
    }

    public int getHeight() {
        return y;
    }

    @Override
    public String toString() {
        return "Board{" +
                "fields=" + Arrays.deepToString(fields) +
                '}';
    }

    public String prettyString() {
        StringBuilder stringBuilder = new StringBuilder();

        for (int y = 0; y < this.y; ++y) {
            for (int x = 0; x < this.x; x++) {
                if (x != 0) stringBuilder.append(" ");
                stringBuilder.append(fields[x][y].prettyString());
            }
            stringBuilder.append("\n");
        }

        return stringBuilder.toString();
    }

    public int countFood() {
        int foodFields = 0;

        for (int y = 0; y < this.y; ++y) {
            for (int x = 0; x < this.x; x++) {
                if (fields[x][y].hasFood()) {
                    ++foodFields;
                }
            }
        }

        return foodFields;
    }

    public void updateFields(int simulationTurn, int renewTime) {
        for (int y = 0; y < this.y; ++y) {
            for (int x = 0; x < this.x; x++) {
                fields[x][y].update(simulationTurn, renewTime);
            }
        }
    }

    public Optional<Direction> whereIsFoodFrom(Coordinates coordinates) {
        for (Direction direction : Direction.values()) {
            Coordinates neighbour = coordinates.from(direction);
            if (getField(neighbour).hasFood()) {
                return Optional.of(direction);
            }
        }

        return Optional.empty();
    }

    public Optional<Coordinates> nearestFood(Coordinates coordinates) {
        List<Coordinates> neighbours = new ArrayList<>(8);

        for (Direction direction : Direction.values()) {
            neighbours.add(coordinates.from(direction));
            neighbours.add(coordinates.from(direction, direction.next()));
        }

        for (Coordinates neighbour : neighbours) {
            if (getField(neighbour).hasFood()) {
                return Optional.of(neighbour);
            }
        }

        return Optional.empty();
    }
}
