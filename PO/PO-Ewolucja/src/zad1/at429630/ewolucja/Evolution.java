package zad1.at429630.ewolucja;

import zad1.at429630.ewolucja.board.Board;
import zad1.at429630.ewolucja.board.Coordinates;
import zad1.at429630.ewolucja.board.Direction;
import zad1.at429630.ewolucja.config.Config;
import zad1.at429630.ewolucja.robs.Rob;
import zad1.at429630.ewolucja.utils.Reporter;

import java.util.ArrayList;
import java.util.List;

public class Evolution {
    private final Config config;
    private final Board board;
    private final List<Rob> robs;
    private final Reporter reporter;

    public Evolution(Config config, Board board) {
        this.config = config;
        this.board = board;
        reporter = new Reporter(board);
        robs = new ArrayList<>();
    }

    public void evolve() {
        System.out.println("Beginning, " + reporter.getReport(robs));

        for (int simulationTurn = 1; simulationTurn <= config.getIleTur(); simulationTurn++) {
            updateBoard(simulationTurn);
            runMoves(simulationTurn);
            reproduce();

            if (isTimeForReport(simulationTurn)) {
                System.out.println(simulationTurn + ", " + reporter.getReport(robs));
            }
        }

        System.out.println("End, " + reporter.getReport(robs));
    }

    private boolean isTimeForReport(int simulationTurn) {
        return simulationTurn % config.getCoIleWypisz() == 0;
    }

    private void updateBoard(int simulationTurn) {
        board.updateFields(simulationTurn, config.getIleRosnieJedzenie());
    }

    public void fillRobs() {
        for (int i = 0; i < config.getPoczIleRobow(); ++i) {
            Coordinates coordinates = Coordinates.randomCoordinates(config);
            robs.add(new Rob(config.getPoczProg(),
                    coordinates,
                    Direction.random(config.getRandomizer()),
                    config.getPoczEnergia(), config, board));
        }
    }

    private void runMoves(int simulationTurn) {
        List<Rob> toBeRemoved = new ArrayList<>();

        for (Rob rob : robs) {
            rob.updateBeforeMove();

            rob.move(coordinates -> board.getField(coordinates).robEntered(rob, simulationTurn));

            rob.updateAfterMove();

            if (!rob.isAlive()) {
                toBeRemoved.add(rob);
            }
        }

        robs.removeAll(toBeRemoved);
    }

    private void reproduce() {
        List<Rob> toBeAdded = new ArrayList<>();

        for (Rob rob : robs) {
            if (rob.canReproduce()) {
                rob.tryToReproduce()
                        .ifPresent(toBeAdded::add);
            }
        }

        robs.addAll(toBeAdded);
    }
}
