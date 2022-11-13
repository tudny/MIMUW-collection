package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

@ToString(callSuper = true)
public class GreaterOrEqual extends NumberComparator {
    public GreaterOrEqual(Instruction argument1, Instruction argument2) {
        super(argument1, argument2);
    }

    @Override
    protected NumberCompareOperation getOperation() {
        return NumberComparator::greaterOrEqual;
    }

    @Override
    protected String getOperationSymbol() {
        return ">=";
    }
}
