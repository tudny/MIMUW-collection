package pl.edu.mimuw.at429630.robson.instructions.loops;

import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.Robson;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Number;
import pl.edu.mimuw.at429630.robson.instructions.Block;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Zero;
import pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic.Minus;
import pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic.Plus;
import pl.edu.mimuw.at429630.robson.instructions.numbers.comparators.GreaterOrEqual;
import pl.edu.mimuw.at429630.robson.instructions.variables.Assignment;
import pl.edu.mimuw.at429630.robson.instructions.variables.Variable;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class WhileTest {
    @Test
    void SimpleWhileTestFromTenToOne() {
        assertDoesNotThrow(() -> {
            RuntimeContext context = new RuntimeContext();

            Instruction begin = new Number(10.0);
            Instruction end = new Number(1.0);

            Instruction counter = new Variable("cnt");
            Instruction acc = new Variable("acc");

            Instruction setup = new Assignment("cnt", begin);
            Instruction setup2 = new Assignment("acc", new Zero());

            Instruction condition = new GreaterOrEqual(counter, end);

            Instruction decrement = new Assignment("cnt", new Minus(new Variable("cnt"), new Number(1.0)));
            Instruction inc = new Assignment("acc", new Plus(new Variable("acc"), new Variable("cnt")));

            Instruction whileBlock = new Block(List.of(inc, decrement));
            Instruction instWhile = new While(condition, whileBlock);

            Instruction block = new Block(List.of(setup, setup2, instWhile, acc));

            Double res = block.execute(context);

            assertEquals(55.0, res);
        });
    }

    @Test
    void TestJsonSum() {
        assertDoesNotThrow(() -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/loops/sum_loop.json");
            Double res = robson.wykonaj();

            assertEquals(55.0, res);
        });
    }
}
