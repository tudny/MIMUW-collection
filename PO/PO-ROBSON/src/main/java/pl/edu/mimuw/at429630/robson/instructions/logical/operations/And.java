package pl.edu.mimuw.at429630.robson.instructions.logical.operations;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.logical.Bool;

@ToString(callSuper = true)
public class And extends LogicOperationTwoArguments {
    public And(Instruction argument1, Instruction argument2) {
        super(argument1, argument2);
    }

    @Override
    protected LogicalOperationFunctionTwoArguments getLogicalOperation() {
        return Bool::and;
    }

    @Override
    protected String getOperationSymbol() {
        return "&&";
    }
}
