package pl.edu.mimuw.at429630.robson.instructions;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Number;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class BlockTest {

    private RuntimeContext context;

    @BeforeEach
    void CreateContext() {
        context = new RuntimeContext();
    }

    @Test
    void EmptyBlockTest() {
        Instruction instruction = Block.EMPTY;

        assertDoesNotThrow(() -> instruction.execute(context));
        Double result = null;
        try {
            result = instruction.execute(context);
        } catch (BladWykonania ignored) { }

        assertEquals(0.0, result);
    }

    @Test
    void ManyNumbersBlock() {
        final int count = 10;
        List<Instruction> instructions = new ArrayList<>(count);

        for (int i = 0; i < count; ++i) {
            instructions.add(new Number((double) i));
        }

        Instruction block = new Block(instructions);

        assertDoesNotThrow(() -> block.execute(context));

        Double result = null;
        try {
            result = block.execute(context);
        } catch (BladWykonania ignored) { }

        assertEquals(((Number) instructions.get(count - 1)).getValue(), result);
    }
}
