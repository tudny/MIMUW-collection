package pl.edu.mimuw.at429630.robson.instructions;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Zero;

import java.util.ArrayList;
import java.util.List;

@ToString
@AllArgsConstructor
public class Block extends Instruction {

    public static final Block EMPTY = new Block();

    @Json(name = "instrukcje")
    private final List<Instruction> instructions;

    private Block() {
        instructions = new ArrayList<>();
    }

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        Double result = 0.0;

        for (Instruction instruction : instructions) {
            result = instruction.execute(context);
        }

        return result;
    }

    @Override
    public String toJava(JavaConverter converter) {
        String name = converter.getNextName(this);

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(name);
        functionBuilder.setToReturn("0.0");

        for (Instruction instruction : instructions) {
            String exeName = instruction.toJava(converter) + "()";
            functionBuilder.addLine(exeName);
            functionBuilder.setToReturn(exeName);
        }

        functionBuilder.removeLastLineIfExists();

        converter.addFunction(functionBuilder);

        return name;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (instructions == null) {
            throwNullArgument();
            return;
        }
        
        for (Instruction instruction : instructions) {
            instruction.validate();
        }
    }
}
