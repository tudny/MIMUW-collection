package pl.edu.mimuw.at429630.robson.instructions.logical.operations;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.TwoArgumentsInstruction;

@ToString
public abstract class LogicOperationTwoArguments extends TwoArgumentsInstruction {

    private static final String CASTER_NAME = "requireBoolean";

    public LogicOperationTwoArguments(Instruction argument1, Instruction argument2) {
        super(argument1, argument2);
    }

    protected abstract LogicalOperationFunctionTwoArguments getLogicalOperation();

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        Double arg1 = argument1.execute(context);
        Double arg2 = argument2.execute(context);

        return getLogicalOperation().getOperation(arg1, arg2);
    }

    @Override
    protected String getCasterName() {
        return CASTER_NAME;
    }
}
