package pl.edu.mimuw.at429630.robson.instructions.loops;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.logical.Bool;

@ToString(callSuper = true)
@AllArgsConstructor
public class While extends Instruction {

    @Json(name = "warunek")
    private final Instruction condition;

    @Json(name = "blok")
    private final Instruction block;

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        while (Bool.javaDoubleToBoolean(condition.execute(context)))
            block.execute(context);

        return 0.0;
    }

    @Override
    public String toJava(JavaConverter converter) {
        String name = converter.getNextName(this);
        String conditionName = condition.toJava(converter);
        String blockName = block.toJava(converter);

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(name);
        functionBuilder.setToReturn("0.0");

        functionBuilder.addLineWithoutSemi(String.format("while (requireBoolean(%s())) {", conditionName));
        functionBuilder.addLine(String.format("\t%s()", blockName));
        functionBuilder.addLineWithoutSemi("}");

        converter.addFunction(functionBuilder);
        
        return name;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (condition == null || block == null) {
            throwNullArgument();
            return;
        }

        condition.validate();
        block.validate();
    }
}
