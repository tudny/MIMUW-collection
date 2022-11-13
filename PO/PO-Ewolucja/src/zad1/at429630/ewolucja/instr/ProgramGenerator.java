package zad1.at429630.ewolucja.instr;

import java.util.ArrayList;
import java.util.List;

public class ProgramGenerator {

    private final List<Instruction> instructions;

    public ProgramGenerator(Program program) {
        this(program.getInstructions());
    }

    public ProgramGenerator(List<Instruction> instructions) {
        this.instructions = new ArrayList<>(instructions);
    }

    public ProgramGenerator() {
        this.instructions = new ArrayList<>();
    }

    public boolean isEmpty() {
        return instructions.isEmpty();
    }

    public int size() {
        return instructions.size();
    }

    public void appendInstruction(Instruction instruction) {
        instructions.add(instruction);
    }

    public void popInstruction() {
        instructions.remove(instructions.size() - 1);
    }

    public void replaceInstruction(Instruction instruction, int index) {
        instructions.remove(index);
        instructions.add(index, instruction);
    }

    public Program makeProgram() {
        return new Program(new ArrayList<>(instructions));
    }
}
