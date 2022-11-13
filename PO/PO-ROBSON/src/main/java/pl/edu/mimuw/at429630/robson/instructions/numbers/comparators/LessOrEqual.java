package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

@ToString(callSuper = true)
public class LessOrEqual extends NumberComparator {
    public LessOrEqual(Instruction argument1, Instruction argument2) {
        super(argument1, argument2);
    }

    @Override
    protected NumberCompareOperation getOperation() {
        return NumberComparator::lessOrEqual;
    }

    @Override
    protected String getOperationSymbol() {
        return "<=";
    }
}
