package pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Number;

import static org.junit.jupiter.api.Assertions.*;

class PlusTest {

    private RuntimeContext context;

    @BeforeEach
    void CreateContext() {
        context = new RuntimeContext();
    }

    Number a, b;
    Instruction instruction;

    @BeforeEach
    void setUp() {
        a = new Number(10.0);
        b = new Number(21.1);

        instruction = new Plus(a, b);
    }

    @Test
    void TestExecute() {
        assertDoesNotThrow(() -> instruction.execute(context));

        try {
            Double result = instruction.execute(context);
            assertEquals(31.1, result, Number.PRECISION);
        } catch (BladWykonania ignored) { }
    }
}
