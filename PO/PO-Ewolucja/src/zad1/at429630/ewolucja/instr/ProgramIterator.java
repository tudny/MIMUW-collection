package zad1.at429630.ewolucja.instr;

public interface ProgramIterator {
    Instruction next();
    boolean hasNext();
    void reset();
}
