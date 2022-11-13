package pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.TwoArgumentsInstruction;

@ToString(callSuper = true)
public abstract class Operation extends TwoArgumentsInstruction {

    public static final String CASTER_NAME = "requireDouble";

    public Operation(Instruction argument1, Instruction argument2) {
        super(argument1, argument2);
    }

    protected abstract NumberOperation getOperation();

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        Double arg1 = argument1.execute(context);
        Double arg2 = argument2.execute(context);

        return getOperation().operation(arg1, arg2);
    }

    @Override
    protected String getCasterName() {
        return CASTER_NAME;
    }
}
