package pl.edu.mimuw.at429630.robson.instructions;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;

@ToString
@AllArgsConstructor
public abstract class TwoArgumentsInstruction extends Instruction {

    @Json(name = "argument1")
    protected final Instruction argument1;

    @Json(name = "argument2")
    protected final Instruction argument2;

    protected abstract String getOperationSymbol();
    protected abstract String getCasterName();

    @Override
    public String toJava(JavaConverter converter) {
        String name = converter.getNextName(this);

        String arg1 = argument1.toJava(converter);
        String arg2 = argument2.toJava(converter);

        String casterName = getCasterName();

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(name);
        functionBuilder.setToReturn(String.format("requireDouble(%s(%s()) %s %s(%s()))",
                casterName, arg1, getOperationSymbol(), casterName, arg2));

        converter.addFunction(functionBuilder);

        return name;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (argument1 == null || argument2 == null) {
            throwNullArgument();
            return;
        }

        argument1.validate();
        argument2.validate();
    }
}
