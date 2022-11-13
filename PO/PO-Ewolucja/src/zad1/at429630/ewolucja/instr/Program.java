package zad1.at429630.ewolucja.instr;

import java.util.List;

public class Program implements ProgramIterator {
    private final List<Instruction> instructions;

    public Program(List<Instruction> instructions) {
        this.instructions = instructions;
    }

    public static Program of(String list) {
        ProgramGenerator programGenerator = new ProgramGenerator();

        for (Instruction instruction : Instruction.arrayOf(list)) {
            programGenerator.appendInstruction(instruction);
        }

        return programGenerator.makeProgram();
    }

    public List<Instruction> getInstructions() {
        return instructions;
    }

    private int currentInstruction = 0;

    @Override
    public Instruction next() {
        int toReturn = currentInstruction;
        currentInstruction = (currentInstruction + 1) % instructions.size();

        return instructions.get(toReturn);
    }

    @Override
    public boolean hasNext() {
        return currentInstruction + 1 < instructions.size();
    }

    @Override
    public void reset() {
        currentInstruction = 0;
    }

    @Override
    public String toString() {
        return "Program{" +
                "instructions=" + instructions +
                '}';
    }
}
