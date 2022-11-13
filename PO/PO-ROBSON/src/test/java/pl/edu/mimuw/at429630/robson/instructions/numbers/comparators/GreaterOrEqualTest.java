package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class GreaterOrEqualTest {
    @Test
    void TestGreaterEqualThan() {
        Double a = 10.0;
        Double b = 20.0;

        assertEquals(false, NumberComparator.greaterOrEqual(a, b));
    }

    @Test
    void TestGreaterEqualThan2() {
        Double a = 10.0;
        Double b = 10.0;

        assertEquals(true, NumberComparator.greaterOrEqual(a, b));
    }

    @Test
    void TestGreaterEqualThan3() {
        Double a = 10.0;
        Double b = 5.0;

        assertEquals(true, NumberComparator.greaterOrEqual(a, b));
    }

    @Test
    void TestGreaterEqualThan4() {
        Double a = -10.0;
        Double b = -20.0;

        assertEquals(true, NumberComparator.greaterOrEqual(a, b));
    }

    @Test
    void TestGreaterEqualThan5() {
        Double a = -10.0;
        Double b = -10.0;

        assertEquals(true, NumberComparator.greaterOrEqual(a, b));
    }

    @Test
    void TestGreaterEqualThan6() {
        Double a = -10.0;
        Double b = -5.0;

        assertEquals(false, NumberComparator.greaterOrEqual(a, b));
    }
}
