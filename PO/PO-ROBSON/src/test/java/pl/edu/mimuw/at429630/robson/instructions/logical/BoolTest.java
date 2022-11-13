package pl.edu.mimuw.at429630.robson.instructions.logical;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.Robson;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;

import static org.junit.jupiter.api.Assertions.*;

class BoolTest {

    private RuntimeContext context;

    @BeforeEach
    void CreateContext() {
        context = new RuntimeContext();
    }

    @Test
    void TestSimpleTrueRead() {
        assertDoesNotThrow(() -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/logical/ok_bool.json");

            double result = robson.wykonaj();

            assertEquals(1.0, result);
        });
    }

    @Test
    void TestWrongDataRead() {
        assertThrows(BladWykonania.class, () -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/logical/wrong_bool.json");
            robson.wykonaj();
        });
    }

    @Test
    void TestBoolStaticConstructor1() {
        assertDoesNotThrow(() -> {
            Bool b = Bool.of(1.0);

            assertEquals(true, b.toBoolean());
        });
    }

    @Test
    void TestBoolStaticConstructor2() {
        assertDoesNotThrow(() -> {
            Bool b = Bool.of(0.0);

            assertEquals(false, b.toBoolean());
        });
    }

    @Test
    void TestBoolStaticConstructorException1() {
        assertThrows(BladWykonania.class, () -> Bool.of(1.1));
    }

    @Test
    void TestBoolStaticConstructorException2() {
        assertThrows(BladWykonania.class, () -> Bool.of(1.00000001));
    }

    @Test
    void TestBoolStaticConstructorException3() {
        assertThrows(BladWykonania.class, () -> Bool.of(0.00000001));
    }

    @Test
    void TestBoolStaticConstructorException4() {
        assertThrows(BladWykonania.class, () -> Bool.of(-0.0));
    }

    @Test
    void TestBoolStaticConstructorNoException() {
        assertDoesNotThrow(() -> Bool.of(-0.0 + 0.0));
    }

    @Test
    void TestExecute1() {
        assertDoesNotThrow(() -> {
            Bool b = Bool.of(1.0);

            Double d = b.execute(context);
            assertEquals(1.0, d);
        });
    }

    @Test
    void TestExecute2() {
        assertDoesNotThrow(() -> {
            Bool b = Bool.of(0.0);

            Double d = b.execute(context);
            assertEquals(0.0, d);
        });
    }

    @Test
    void TestExecuteException() {
        assertThrows(BladWykonania.class, () -> Bool.of(1.1));
    }
}
