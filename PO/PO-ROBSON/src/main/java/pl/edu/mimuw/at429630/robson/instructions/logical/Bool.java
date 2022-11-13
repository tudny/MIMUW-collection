package pl.edu.mimuw.at429630.robson.instructions.logical;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

import java.util.Objects;

@ToString
public abstract class Bool extends Instruction {
    private final Double value;

    protected Bool(Double value) {
        this.value = value;
    }

    public static Bool of(Double x) throws BladWykonania {
        return javaDoubleToBoolean(x) ? new True() : new False();
    }

    public static Boolean javaDoubleToBoolean(Double value) throws BladWykonania {
        if (value.equals(0.0))
            return false;
        else if (value.equals(1.0))
            return true;

        throw new BladWykonania("Trying to convert " + value + " to Boolean and this is not bool bro.");
    }

    public static Double javaBooleanToDouble(Boolean value) {
        return value ? 1.0 : 0.0;
    }

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        return value;
    }

    public static Double operationTwoArguments(Double a, Double b, BoolOperationTwoArguments boolOperationTwoArguments)
            throws BladWykonania {
        Boolean arg1 = javaDoubleToBoolean(a);
        Boolean arg2 = javaDoubleToBoolean(b);

        return javaBooleanToDouble(boolOperationTwoArguments.operation(arg1, arg2));
    }

    public static Double operationSingleArgument(Double a, BoolOperationSingleArgument boolOperationSingleArgument)
            throws BladWykonania {
        Boolean arg = javaDoubleToBoolean(a);

        return javaBooleanToDouble(boolOperationSingleArgument.operation(arg));
    }

    public static Double and(Double a, Double b) throws BladWykonania {
        return operationTwoArguments(a, b, (a1, b1) -> a1 && b1);
    }

    public static Double or(Double a, Double b) throws BladWykonania {
        return operationTwoArguments(a, b, (a1, b1) -> a1 || b1);
    }

    public static Double not(Double a) throws BladWykonania {
        return operationSingleArgument(a, (a1) -> !a1);
    }

    public Boolean toBoolean() throws BladWykonania {
        return javaDoubleToBoolean(value);
    }

    @FunctionalInterface
    private interface BoolOperationTwoArguments {
        Boolean operation(Boolean a, Boolean b);
    }

    @FunctionalInterface
    private interface BoolOperationSingleArgument {
        Boolean operation(Boolean a);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null) return false;
        if (!(o instanceof Bool)) return false;
        Bool bool = (Bool) o;
        return Objects.equals(value, bool.value);
    }

    @Override
    public String toJava(JavaConverter converter) {
        String name = converter.getNextName(this);

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(name);
        functionBuilder.setToReturn(value.toString());

        converter.addFunction(functionBuilder);

        return name;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (value == null)
            throwNullArgument();
    }
}
