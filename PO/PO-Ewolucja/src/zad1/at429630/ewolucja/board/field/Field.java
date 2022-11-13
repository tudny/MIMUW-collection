package zad1.at429630.ewolucja.board.field;

import zad1.at429630.ewolucja.robs.Rob;

public abstract class Field {

    public static Field of(Character character) {
        switch (character) {
            case ' ':
                return new EmptyField();
            case 'x':
                return new FoodField();
            default:
                throw new IllegalStateException("Unexpected value: '" + character + "' during Field creation.");
        }
    }

    public void robEntered(Rob rob, int actualTime) { }

    public abstract void update(int actualTime, int renewTime);

    public abstract boolean hasFood();

    @Override
    public String toString() {
        return "Field{}";
    }

    public abstract String prettyString();
}
