package pl.edu.mimuw.at429630.robson.instructions.numbers;

import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.Robson;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;

import static org.junit.jupiter.api.Assertions.*;

class NumberTest {
    @Test
    void TestSimpleNumberRead() {
        assertDoesNotThrow(() -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/number/simple_number.json");

            double result = robson.wykonaj();

            assertEquals(42, result);
        });
    }

    @Test
    void TestWrongNumberRead() {
        Robson robson = new Robson();
        assertThrows(NieprawidlowyProgram.class, () -> robson.fromJSON("tests/number/wrong_number.json"));
    }

    @Test
    void TestSum1() {
        Number a = new Number(10.0);
        Number b = new Number(20.0);

        Double expected = 30.0;
        Double result = Number.sum(a.getValue(), b.getValue());

        assertEquals(expected, result);
    }

    @Test
    void TestSum2() {
        Number a = new Number((double) (1 << 20));
        Number b = new Number((double) -(1 << 20));

        Double expected = 0.0;
        Double result = Number.sum(a.getValue(), b.getValue());

        assertEquals(expected, result);
    }

    @Test
    void TestSumLong() {
        Number base = new Zero();

        for (double a = 0.0; a <= 10.0; a += 1.0) {
            Number aNumber = new Number(a);
            base = new Number(Number.sum(aNumber.getValue(), base.getValue()));
        }

        double start = 0.0;
        double end = 10.0;
        Number expected = new Number((start + end) * (end - start + 1.0) / 2.0);

        assertEquals(expected, base);
    }

    @Test
    void TestDifference1() {
        Number a = new Number(42.0);
        Number b = new Number(720.0);

        Double expected = -678.0;
        Double result = Number.difference(a.getValue(), b.getValue());

        assertEquals(expected, result);
    }

    @Test
    void TestDifference2() {
        Number a = new Number(-123.0);
        Number b = new Number(-123.0);

        Double expected = 0.0;
        Double result = Number.difference(a.getValue(), b.getValue());

        assertEquals(expected, result);
    }

    @Test
    void TestProduct1() {
        Number a = new Number(-123.0);
        Number b = new Number(-123.0);

        Double expected = 15129.0;
        Double result = Number.product(a.getValue(), b.getValue());

        assertEquals(expected, result);
    }

    @Test
    void TestProduct2() {
        Number a = new Number(0.0);
        Number b = new Number(-123.0);

        Double expected = 0.0;
        Double result = Number.product(a.getValue(), b.getValue());

        assertEquals(expected, result);
    }

    @Test
    void TestProductLong() {
        Number base = new One();
        double expected = 1.0;

        for (double a = 1.0; a <= 10.0; a += 1.0) {
            Number aNumber = new Number(a);
            base = new Number(Number.product(base.getValue(), aNumber.getValue()));

            expected *= a;
        }

        Number expectedNumber = new Number(expected);

        assertEquals(expectedNumber, base);
    }

    @Test
    void TestQuotient1() {
        Number a = new Number(20.0);
        Number b = new Number(2.0);

        Double expected = 10.0;
        Double result = Number.quotient(a.getValue(), b.getValue());

        assertEquals(expected, result);
    }

    @Test
    void TestQuotient2() {
        Number a = new Number(0.0);
        Number b = new Number(2.0);

        Double expected = 0.0;
        Double result = Number.quotient(a.getValue(), b.getValue());

        assertEquals(expected, result);
    }

    @Test
    void TestQuotientZeroDivision() {
        Number a = new Number(1.0);
        Number b = new Zero();

        assertThrows(ArithmeticException.class, () -> Number.quotient(a.getValue(), b.getValue()));
    }

    @Test
    void TestIsZeroConstructor() {
        Number a = new Number(0.0);
        Number b = new Zero();

        assertEquals(a, b);
        assertEquals(b, a);

        assertEquals(Boolean.TRUE, a.isZero());
        assertEquals(Boolean.TRUE, b.isZero());
    }

    @Test
    void TestIsZeroAfterOperation1() {
        Number a = new Number(2137.0);
        Number b = new Zero();

        Double result = Number.product(a.getValue(), b.getValue());

        assertEquals(0.0, result);
    }

    @Test
    void TestIsZeroAfterOperation2() {
        Number a = new Number(-2137.0);
        Number b = new Zero();

        Double result = Number.product(a.getValue(), b.getValue());

        assertEquals(0.0, result);
    }

    @Test
    void TestIsOneConstructor() {
        Number a = new Number(1.0);
        Number b = new One();

        assertEquals(a, b);
        assertEquals(b, a);

        assertEquals(Boolean.TRUE, a.isOne());
        assertEquals(Boolean.TRUE, b.isOne());
    }

    @Test
    void TestIsOneAfterOperation() {
        Number a = new Number(2137.0);
        Number b = new Number(2137.0);

        Double result = Number.quotient(a.getValue(), b.getValue());
        Double expected = 1.0;

        assertEquals(expected, result);
    }
}
