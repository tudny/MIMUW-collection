package pl.edu.mimuw.at429630.robson.instructions.logical.operations;

import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.Robson;

import static org.junit.jupiter.api.Assertions.*;

class AndTest {
    @Test
    void TestSimpleOperationRead() {
        assertDoesNotThrow(() -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/logical/operations/simple_and.json");

            double result = robson.wykonaj();

            assertEquals(0.0, result);
        });
    }

    @Test
    void TestSimpleOperationRead2() {
        assertDoesNotThrow(() -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/logical/operations/simple_and_2.json");

            double result = robson.wykonaj();

            assertEquals(1.0, result);
        });
    }

    @Test
    void TestSimpleOperationRead3() {
        assertDoesNotThrow(() -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/logical/operations/simple_and_3.json");

            double result = robson.wykonaj();

            assertEquals(0.0, result);
        });
    }

    @Test
    void TestComplexOperationRead() {
        assertDoesNotThrow(() -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/logical/operations/complex_and.json");

            double result = robson.wykonaj();

            assertEquals(1.0, result);
        });
    }
}
