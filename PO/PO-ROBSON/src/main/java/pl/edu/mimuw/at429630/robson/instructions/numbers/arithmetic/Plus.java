package pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Number;

@ToString(callSuper = true)
public class Plus extends Operation {

    public Plus(Instruction argument1, Instruction argument2) {
        super(argument1, argument2);
    }

    @Override
    protected NumberOperation getOperation() {
        return Number::sum;
    }

    @Override
    protected String getOperationSymbol() {
        return "+";
    }
}
