package pl.edu.mimuw.at429630.robson.instructions.logical.operations;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.logical.Bool;

@ToString(callSuper = true)
public class Not extends LogicOperationSingleArgument {

    public Not(Instruction argument) {
        super(argument);
    }

    @Override
    protected LogicalOperationFunctionSingleArgument getLogicalOperation() {
        return Bool::not;
    }

    @Override
    protected String getOperationSymbol() {
        return "!";
    }
}
