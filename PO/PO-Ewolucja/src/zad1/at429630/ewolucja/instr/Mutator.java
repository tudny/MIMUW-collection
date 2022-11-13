package zad1.at429630.ewolucja.instr;

import zad1.at429630.ewolucja.config.Config;

public class Mutator {
    private final Config config;

    public Mutator(Config config) {
        this.config = config;
    }

    public Program mutate(Program program) {
        ProgramGenerator programGenerator = new ProgramGenerator(program);

        tryRemovingLast(programGenerator);
        tryAddingNext(programGenerator);
        tryChangingOne(programGenerator);

        return programGenerator.makeProgram();
    }

    private Instruction randomInstruction() {
        Instruction[] availableInstructions = config.getSpisInstr();
        int rng = config.getRandomizer().getRng(availableInstructions.length);
        return availableInstructions[rng];
    }

    private void tryRemovingLast(ProgramGenerator programGenerator) {
        if (config.getRandomizer().happens(config.getPrUsunieciaInstr())) {
            if (!programGenerator.isEmpty()) {
                programGenerator.popInstruction();
            }
        }
    }

    private void tryAddingNext(ProgramGenerator programGenerator) {
        if (config.getRandomizer().happens(config.getPrDodanieInstr())) {
            programGenerator.appendInstruction(randomInstruction());
        }
    }

    private void tryChangingOne(ProgramGenerator programGenerator) {
        if (config.getRandomizer().happens(config.getPrZmianaInstr())) {
            if (!programGenerator.isEmpty()) {
                int randomIndex = config.getRandomizer().getRng(programGenerator.size());
                programGenerator.replaceInstruction(randomInstruction(), randomIndex);
            }
        }
    }
}
