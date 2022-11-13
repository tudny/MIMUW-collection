package zad1.at429630.ewolucja.robs;

import zad1.at429630.ewolucja.board.Board;
import zad1.at429630.ewolucja.board.Coordinates;
import zad1.at429630.ewolucja.board.Direction;
import zad1.at429630.ewolucja.config.Config;
import zad1.at429630.ewolucja.instr.Instruction;
import zad1.at429630.ewolucja.instr.Mutator;
import zad1.at429630.ewolucja.instr.Program;
import zad1.at429630.ewolucja.instr.ProgramIterator;

import java.util.Optional;

public class Rob {
    private final static int INSTRUCTION_EXECUTION_COST = 1;

    private final Program program;
    private final Board board;
    private final Config config;
    private Coordinates coordinates;
    private Direction rotation;
    private int energy;
    private int age;

    public Rob(Program program, Coordinates coordinates, Direction direction, int beginEnergy, Config config, Board board) {
        this.program = program;
        this.coordinates = coordinates;
        this.energy = beginEnergy;
        this.config = config;
        this.rotation = direction;
        this.board = board;
        this.age = 0;
    }

    public int getProgramLength() {
        return program.getInstructions().size();
    }

    public int getEnergy() {
        return energy;
    }

    public int getAge() {
        return age;
    }

    public void updateBeforeMove() {
        ++age;
    }

    public void move(OnCoordinatesChangeListener onCoordinatesChangeListener) {
        ProgramIterator programIterator = program;
        while (programIterator.hasNext()) {
            if (!canMove()) {
                break;
            }

            Instruction instruction = programIterator.next();
            executeInstruction(instruction, onCoordinatesChangeListener);
        }

        programIterator.reset();
    }

    public void updateAfterMove() {
        takeEnergy(config.getKosztTury());
    }

    public boolean isAlive() {
        return energy >= 0;
    }

    public boolean canMove() {
        return energy >= 1;
    }

    public boolean canReproduce() {
        return energy >= config.getLimitPowielania();
    }

    public void takeEnergy(int value) {
        energy -= value;
    }

    public void addEnergy(int value) {
        energy += value;
    }

    public void addDefaultEnergy() {
        addEnergy(config.getIleDajeJedzenia());
    }

    private void executeInstruction(Instruction instruction,
                                    OnCoordinatesChangeListener onCoordinatesChangeListener) {
        takeEnergy(INSTRUCTION_EXECUTION_COST);

        switch (instruction) {
            case LEFT:
                rotation = rotation.next();
                break;
            case RIGHT:
                rotation = rotation.prev();
                break;
            case FORWARD:
                moveTo(coordinates.from(rotation), onCoordinatesChangeListener);
                break;
            case SMELL:
                findFood();
                break;
            case EAT:
                goForFood(onCoordinatesChangeListener);
                break;
        }
    }

    private void findFood() {
        board.whereIsFoodFrom(coordinates)
                .ifPresent(direction -> this.rotation = direction);
    }

    private void goForFood(OnCoordinatesChangeListener onCoordinatesChangeListener) {
        board.nearestFood(coordinates)
                .ifPresent(coordinates -> moveTo(coordinates, onCoordinatesChangeListener));
    }

    private void moveTo(Coordinates coordinates, OnCoordinatesChangeListener onCoordinatesChangeListener) {
        this.coordinates = coordinates;
        onCoordinatesChangeListener.coordinatesChangedCallback(coordinates);
    }

    public Optional<Rob> tryToReproduce() {
        if (config.getRandomizer().happens(config.getPrPowielania())) {
            Program newProgram = new Mutator(config).mutate(program);
            int newEnergy = (int) (energy * config.getUlamekEnergiiRodzica());
            Coordinates newCoordinates = coordinates;
            Direction newDirection = rotation.next().next();

            takeEnergy(newEnergy);
            return Optional.of(new Rob(newProgram, newCoordinates, newDirection, newEnergy, config, board));
        }

        return Optional.empty();
    }

    @Override
    public String toString() {
        return "Rob{" +
                "program=" + program +
                ", coordinates=" + coordinates +
                ", rotation=" + rotation +
                ", energy=" + energy +
                '}';
    }
}
