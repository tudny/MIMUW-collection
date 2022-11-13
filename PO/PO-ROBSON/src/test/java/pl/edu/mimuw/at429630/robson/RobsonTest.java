package pl.edu.mimuw.at429630.robson;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class RobsonTest {
    private Robson robson;

    @BeforeEach
    void Setup() {
        robson = new Robson();
    }

    @Test
    void TestOfficialSum() {
        assertDoesNotThrow(() -> {
            robson.fromJSON("tests/official/add.json");
            Double res = robson.wykonaj();

            assertEquals(15.0, res);
        });
    }

    @Test
    void TestOfficialFib() {
        assertDoesNotThrow(() -> {
            robson.fromJSON("tests/official/fibb.json");
            Double res = robson.wykonaj();

            assertEquals(55.0, res);
        });
    }

    @Test
    void TestToJavaSum() {
        assertDoesNotThrow(() -> {
            robson.fromJSON("tests/official/add.json");
            robson.toJava("add.java");

            System.out.println("The result should be checked by the user. file add.java");
        });
    }

    @Test
    void TestAllInOne() {
        assertDoesNotThrow(() -> {
            robson.fromJSON("tests/logical/not_java.json");
            robson.toJava("not.java");
            Double res = robson.wykonaj();
            assertEquals((double) (8 + 6 + 4 + 2), res);

            System.out.println("The result should be checked by the user. file not.java");
        });
    }
}
