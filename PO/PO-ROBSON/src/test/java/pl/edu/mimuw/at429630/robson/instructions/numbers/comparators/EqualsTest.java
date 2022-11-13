package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class EqualsTest {
    @Test
    void TestEqual() {
        Double a = 1.0;
        Double b = 1.0;

        assertTrue(NumberComparator.equals(a, b));
    }

    @Test
    void TestEqual2() {
        Double a = 1.0;
        Double b = 0.1;

        assertFalse(NumberComparator.equals(a, b));
    }

    @Test
    void TestEqualCornerCase() {
        Double a = -0.0;
        Double b = 0.0;

        assertFalse(NumberComparator.equals(a, b));
    }

    @Test
    void TestEqualCornerCase2() {
        Double a = Double.NaN;
        Double b = 1.123123123123123;

        assertFalse(NumberComparator.equals(a, b));
    }

    @Test
    void TestEqualCornerCase3() {
        Double a = Double.MAX_VALUE;
        Double b = Double.MIN_VALUE;

        assertFalse(NumberComparator.equals(a, b));
    }
}
