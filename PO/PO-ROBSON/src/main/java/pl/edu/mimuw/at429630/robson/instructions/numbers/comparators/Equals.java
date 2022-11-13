package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Number;

import java.util.Locale;

@ToString(callSuper = true)
public class Equals extends NumberComparator {
    public Equals(Instruction argument1, Instruction argument2) {
        super(argument1, argument2);
    }

    @Override
    protected NumberCompareOperation getOperation() {
        return NumberComparator::equals;
    }

    @Override
    protected String getOperationSymbol() {
        return "==";
    }

    @Override
    public String toJava(JavaConverter converter) {
        String name = converter.getNextName(this);

        String arg1 = argument1.toJava(converter);
        String arg2 = argument2.toJava(converter);

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(name);
        functionBuilder.setToReturn(String.format(Locale.US, "requireDouble(Math.abs(%s() - %s()) < %.10f)",
                arg1, arg2, Number.PRECISION));

        converter.addFunction(functionBuilder);

        return name;
    }
}
