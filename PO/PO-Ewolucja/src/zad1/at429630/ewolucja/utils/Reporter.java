package zad1.at429630.ewolucja.utils;

import zad1.at429630.ewolucja.board.Board;
import zad1.at429630.ewolucja.robs.Rob;

import java.util.List;

public class Reporter {
    private final Board board;

    public Reporter(Board board) {
        this.board = board;
    }

    public String getReport(List<Rob> robs) {
        int food = board.countFood();

        String reportString = String.format("rob: %d, żyw: %d ", robs.size(), food);

        String program = String.format("prg: %s, ", reportBy(Rob::getProgramLength, robs));
        String energy = String.format("energ: %s, ", reportBy(Rob::getEnergy, robs));
        String age = String.format("wiek: %s", reportBy(Rob::getAge, robs));

        return reportString + program + energy + age;
    }

    private String reportBy(DataCollector dataCollector, List<Rob> robs) {
        if (robs.isEmpty()) {
            return "0/0/0";
        }

        int minimal = Integer.MAX_VALUE;
        int sum = 0;
        int maximal = Integer.MIN_VALUE;
        int count = robs.size();

        for (Rob rob : robs) {
            minimal = Math.min(minimal, dataCollector.getData(rob));
            sum += dataCollector.getData(rob);
            maximal = Math.max(maximal, dataCollector.getData(rob));
        }

        double average = sum / (double) count;

        return String.format("%d/%.2f/%d", minimal, average, maximal);
    }
}
