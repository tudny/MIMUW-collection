package pl.edu.mimuw.at429630.robson.instructions.logical.operations;

import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.Robson;

import static org.junit.jupiter.api.Assertions.*;

class NotTest {
    @Test
    void TestSimpleOperationRead() {
        assertDoesNotThrow(() -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/logical/operations/simple_not.json");

            double result = robson.wykonaj();

            assertEquals(0.0, result);
        });
    }
}
