package pl.edu.mimuw.at429630.robson.instructions;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;

@ToString
@AllArgsConstructor
public abstract class SingleArgumentPreInstruction extends Instruction {

    @Json(name = "argument")
    protected final Instruction argument;

    protected abstract String getOperationSymbol();
    protected abstract String getCasterName();

    @Override
    public String toJava(JavaConverter converter) {
        String name = converter.getNextName(this);
        String arg = argument.toJava(converter);
        String casterName = getCasterName();

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(name);
        functionBuilder.setToReturn(String.format("requireDouble(%s%s(%s()))", getOperationSymbol(), casterName, arg));

        converter.addFunction(functionBuilder);

        return name;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (argument == null) {
            throwNullArgument();
            return;
        }

        argument.validate();
    }
}
