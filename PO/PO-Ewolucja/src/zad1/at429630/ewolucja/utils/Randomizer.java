package zad1.at429630.ewolucja.utils;

import java.util.Random;

public class Randomizer {
    private final Random rand = new Random();

    public Randomizer() {
        rand.setSeed(System.currentTimeMillis());
    }

    public boolean happens(double chance) {
        return rand.nextDouble() < chance;
    }

    public int getRng(int upperBound) {
        return rand.nextInt(upperBound);
    }
}
