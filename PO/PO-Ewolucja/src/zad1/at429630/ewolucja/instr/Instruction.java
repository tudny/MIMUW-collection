package zad1.at429630.ewolucja.instr;

import java.util.HashSet;
import java.util.Set;

public enum Instruction {
    LEFT('l'),
    RIGHT('p'),
    FORWARD('i'),
    SMELL('w'),
    EAT('j');

    private final Character id;

    Instruction(Character id) {
        this.id = id;
    }

    public static Instruction of(Character character) {
        for (Instruction instruction : Instruction.values()) {
            if (instruction.id.equals(character))
                return instruction;
        }

        return null;
    }

    public static Instruction[] arrayOf(String str) {
        Set<Instruction> instructions = new HashSet<>();

        for (Character character : str.toCharArray()) {
            instructions.add(Instruction.of(character));
        }

        return instructions.toArray(new Instruction[0]);
    }
}
