package pl.edu.mimuw.at429630.robson.instructions.variables;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Number;
import pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic.Plus;

import static org.junit.jupiter.api.Assertions.*;

class VariableAssignmentTest {

    private RuntimeContext context;

    @BeforeEach
    void SetupContext() {
        context = new RuntimeContext();
    }

    @Test
    void TestSimpleAssignment() {
        assertDoesNotThrow(() -> {
            Instruction number = new Number(10.0);
            Instruction instruction = new Assignment("abc", number);

            Double afterAssign = instruction.execute(context);

            assertEquals(10.0, afterAssign);
        });
    }

    @Test
    void TestComplexAssignment() {
        assertDoesNotThrow(() -> {
            Instruction number1 = new Number(10.0);
            Instruction number2 = new Number(123.123);
            Instruction sum = new Plus(number1, number2);

            Instruction instruction = new Assignment("abc", sum);

            Double afterAssign = instruction.execute(context);

            assertEquals(133.123, afterAssign);
        });
    }

    @Test
    void TestDefaultVariableValue() {
        assertDoesNotThrow(() -> {
            Instruction instruction = new Variable("abc");
            Double res = instruction.execute(context);

            assertEquals(0.0, res);
        });
    }
}
