package zad1.at429630.ewolucja.board.field;

import zad1.at429630.ewolucja.robs.Rob;

public class FoodField extends Field {

    private int lastFoodEaten = 0;
    private boolean hasFood = true;

    @Override
    public void robEntered(Rob rob, int actualTime) {
        super.robEntered(rob, actualTime);

        if (hasFood) {
            lastFoodEaten = actualTime;
            hasFood = false;
            rob.addDefaultEnergy();
        }
    }

    @Override
    public void update(int actualTime, int renewTime) {
        if (actualTime == lastFoodEaten + renewTime) {
            hasFood = true;
        }
    }

    @Override
    public boolean hasFood() {
        return hasFood;
    }

    @Override
    public String toString() {
        return "FoodField{" +
                "lastFoodEaten=" + lastFoodEaten +
                ", hasFood=" + hasFood +
                '}';
    }

    @Override
    public String prettyString() {
        return hasFood ? "X" : "x";
    }
}
