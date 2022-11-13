package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.TwoArgumentsInstruction;
import pl.edu.mimuw.at429630.robson.instructions.logical.Bool;

@ToString(callSuper = true)
public abstract class NumberComparator extends TwoArgumentsInstruction {

    private static final String CASTER_NAME = "requireDouble";

    public NumberComparator(Instruction argument1, Instruction argument2) {
        super(argument1, argument2);
    }

    protected abstract NumberCompareOperation getOperation();

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        Double arg1 = argument1.execute(context);
        Double arg2 = argument2.execute(context);

        return Bool.javaBooleanToDouble(getOperation().operation(arg1, arg2));
    }

    public static Boolean equals(Double a, Double b) {
        return a.equals(b);
    }

    public static Boolean lessThan(Double a, Double b) {
        return a < b;
    }

    public static Boolean greaterThan(Double a, Double b) {
        return a > b;
    }

    public static Boolean lessOrEqual(Double a, Double b) {
        return lessThan(a,b) || equals(a, b);
    }

    public static Boolean greaterOrEqual(Double a, Double b) {
        return greaterThan(a,b) || equals(a, b);
    }

    @Override
    protected String getCasterName() {
        return CASTER_NAME;
    }
}
