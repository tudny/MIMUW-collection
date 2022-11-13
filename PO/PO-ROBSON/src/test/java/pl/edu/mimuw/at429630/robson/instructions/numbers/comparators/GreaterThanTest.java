package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class GreaterThanTest {
    @Test
    void TestGreaterThan() {
        Double a = 10.0;
        Double b = 20.0;

        assertEquals(false, NumberComparator.greaterThan(a, b));
    }

    @Test
    void TestGreaterThan2() {
        Double a = 10.0;
        Double b = 10.0;

        assertEquals(false, NumberComparator.greaterThan(a, b));
    }

    @Test
    void TestGreaterThan3() {
        Double a = 10.0;
        Double b = 5.0;

        assertEquals(true, NumberComparator.greaterThan(a, b));
    }

    @Test
    void TestGreaterThan4() {
        Double a = -10.0;
        Double b = -20.0;

        assertEquals(true, NumberComparator.greaterThan(a, b));
    }

    @Test
    void TestGreaterThan5() {
        Double a = -10.0;
        Double b = -10.0;

        assertEquals(false, NumberComparator.greaterThan(a, b));
    }

    @Test
    void TestGreaterThan6() {
        Double a = -10.0;
        Double b = -5.0;

        assertEquals(false, NumberComparator.greaterThan(a, b));
    }
}
