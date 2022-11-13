package pl.edu.mimuw.at429630.robson.instructions.logical;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

@ToString
@AllArgsConstructor
public class If extends Instruction {

    @Json(name = "warunek")
    private final Instruction condition;

    @Json(name = "blok_prawda")
    private final Instruction trueBlock;

    @Json(name = "blok_falsz")
    private final Instruction falseBlock;

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        Double cond = condition.execute(context);

        if (Bool.javaDoubleToBoolean(cond))
            return trueBlock.execute(context);

        return falseBlock != null ? falseBlock.execute(context) : 0.0;
    }

    @Override
    public String toJava(JavaConverter converter) {
        String name = converter.getNextName(this);

        String conditionStr = condition.toJava(converter);
        String trueStr = trueBlock.toJava(converter);
        String falseStr = "";
        if (falseBlock != null)
            falseStr = falseBlock.toJava(converter);

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(name);

        functionBuilder.addLine("Double result = 0.0");
        functionBuilder.addLineWithoutSemi(String.format("if (requireBoolean(%s())) {", conditionStr));
        functionBuilder.addLine(String.format("\tresult = %s()", trueStr));

        if (falseBlock != null) {
            functionBuilder.addLineWithoutSemi("} else {");
            functionBuilder.addLine(String.format("\tresult = %s()", falseStr));
        }

        functionBuilder.addLineWithoutSemi("}");
        functionBuilder.setToReturn("result");

        converter.addFunction(functionBuilder);

        return name;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (trueBlock == null || condition == null) {
            throwNullArgument();
            return;
        }

        trueBlock.validate();
        condition.validate();
        if (falseBlock != null)
            falseBlock.validate();
    }
}
