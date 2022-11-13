package pl.edu.mimuw.at429630.robson.instructions.variables;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

@ToString(callSuper = true)
@AllArgsConstructor
public class Variable extends Instruction {

    @Json(name = "nazwa")
    private final String name;

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        return context.getVariable(name);
    }

    @Override
    public String toJava(JavaConverter converter) {
        String functionName = converter.getNextName(this);

        FunctionBuilder functionBuilder = new FunctionBuilder();

        functionBuilder.setName(functionName);
        functionBuilder.setToReturn(this.name);

        converter.addFunction(functionBuilder);

        return functionName;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (name == null || name.equals(""))
            throwNullArgument();
    }
}
