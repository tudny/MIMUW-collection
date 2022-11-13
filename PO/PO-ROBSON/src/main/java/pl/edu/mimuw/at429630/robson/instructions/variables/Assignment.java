package pl.edu.mimuw.at429630.robson.instructions.variables;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.conversion.VariableBuilder;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

@ToString(callSuper = true)
@AllArgsConstructor
public class Assignment extends Instruction {

    @Json(name = "nazwa")
    private final String name;

    @Json(name = "wartosc")
    private final Instruction value;

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        Double toAssign = value.execute(context);

        context.setVariable(name, toAssign);

        return toAssign;
    }

    @Override
    public String toJava(JavaConverter converter) {
        String functionName = converter.getNextName(this);
        String toAssignName = value.toJava(converter);

        converter.addVariable(new VariableBuilder().name(name));

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(functionName);
        functionBuilder.setToReturn(String.format("%s = %s()", name, toAssignName));

        converter.addFunction(functionBuilder);

        return functionName;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (name == null || name.equals("") || value == null) {
            throwNullArgument();
            return;
        }

        value.validate();
    }
}
