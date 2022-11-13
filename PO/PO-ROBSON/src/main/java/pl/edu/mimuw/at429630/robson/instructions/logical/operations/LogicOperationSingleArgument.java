package pl.edu.mimuw.at429630.robson.instructions.logical.operations;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.SingleArgumentPreInstruction;

@ToString
public abstract class LogicOperationSingleArgument extends SingleArgumentPreInstruction {
    private static final String CASTER_NAME = "requireBoolean";

    public LogicOperationSingleArgument(Instruction argument) {
        super(argument);
    }

    protected abstract LogicalOperationFunctionSingleArgument getLogicalOperation();

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        Double arg = argument.execute(context);

        return getLogicalOperation().getOperation(arg);
    }

    @Override
    protected String getCasterName() {
        return CASTER_NAME;
    }
}
