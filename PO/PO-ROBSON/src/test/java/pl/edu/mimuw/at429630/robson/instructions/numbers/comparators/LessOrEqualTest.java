package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class LessOrEqualTest {
    @Test
    void TestLessEqualThan() {
        Double a = 10.0;
        Double b = 20.0;

        assertEquals(true, NumberComparator.lessOrEqual(a, b));
    }

    @Test
    void TestLessEqualThan2() {
        Double a = 10.0;
        Double b = 10.0;

        assertEquals(true, NumberComparator.lessOrEqual(a, b));
    }

    @Test
    void TestLessEqualThan3() {
        Double a = 10.0;
        Double b = 5.0;

        assertEquals(false, NumberComparator.lessOrEqual(a, b));
    }

    @Test
    void TestLessEqualThan4() {
        Double a = -10.0;
        Double b = -20.0;

        assertEquals(false, NumberComparator.lessOrEqual(a, b));
    }

    @Test
    void TestLessEqualThan5() {
        Double a = -10.0;
        Double b = -10.0;

        assertEquals(true, NumberComparator.lessOrEqual(a, b));
    }

    @Test
    void TestLessEqualThan6() {
        Double a = -10.0;
        Double b = -5.0;

        assertEquals(true, NumberComparator.lessOrEqual(a, b));
    }
}
