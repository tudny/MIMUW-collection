package zad1.at429630.ewolucja.board.field;

public class EmptyField extends Field {

    @Override
    public void update(int actualTime, int renewTime) { }

    @Override
    public boolean hasFood() {
        return false;
    }

    @Override
    public String toString() {
        return "EmptyField{}";
    }

    @Override
    public String prettyString() {
        return ".";
    }
}
