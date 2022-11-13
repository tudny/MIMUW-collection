package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class NumberComparatorTest {
    @Test
    void TestLessThan() {
        Double a = 10.0;
        Double b = 20.0;

        assertEquals(true, NumberComparator.lessThan(a, b));
    }

    @Test
    void TestLessThan2() {
        Double a = 10.0;
        Double b = 10.0;

        assertEquals(false, NumberComparator.lessThan(a, b));
    }

    @Test
    void TestLessThan3() {
        Double a = 10.0;
        Double b = 5.0;

        assertEquals(false, NumberComparator.lessThan(a, b));
    }

    @Test
    void TestLessThan4() {
        Double a = -10.0;
        Double b = -20.0;

        assertEquals(false, NumberComparator.lessThan(a, b));
    }

    @Test
    void TestLessThan5() {
        Double a = -10.0;
        Double b = -10.0;

        assertEquals(false, NumberComparator.lessThan(a, b));
    }

    @Test
    void TestLessThan6() {
        Double a = -10.0;
        Double b = -5.0;

        assertEquals(true, NumberComparator.lessThan(a, b));
    }
}
