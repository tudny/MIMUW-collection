package pl.edu.mimuw.at429630.robson.instructions.numbers;

import com.squareup.moshi.Json;
import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.FunctionBuilder;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

@ToString
@AllArgsConstructor
public class Number extends Instruction {

    public static final Double PRECISION = 1e-10;
    public static final Double DEFAULT_GLOBAL_VALUE = 0.0;

    @Json(name = "wartosc")
    @Getter
    private final Double value;

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        return value;
    }

    // Dodawanie 0.0 do każdego wyniku ma na calu zniwelowanie rozwiązań -0.0 do 0.0

    public static Double sum(Double a, Double b) {
        return a + b + 0.0;
    }

    public static Double difference(Double a, Double b) {
        return a - b + 0.0;
    }

    public static Double product(Double a, Double b) {
        return a * b + 0.0;
    }

    public static Double quotient(Double a, Double b) {
        if (b.equals(0.0)) throw new ArithmeticException("Dividing by zero is not allowed!");
        return a / b + 0.0;
    }

    public boolean isZero() {
        return this.equals(new Zero());
    }

    public boolean isOne() {
        return this.equals(new One());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null) return false;
        if (!(o instanceof Number)) return false;
        Number number = (Number) o;
        return Math.abs(value - number.value) < PRECISION;
    }

    @Override
    public String toJava(JavaConverter converter) {
        String name = converter.getNextName(this);

        FunctionBuilder functionBuilder = new FunctionBuilder();
        functionBuilder.setName(name);
        functionBuilder.setToReturn(value.toString());

        converter.addFunction(functionBuilder);

        return name;
    }

    @Override
    public void validate() throws NieprawidlowyProgram {
        if (value == null)
            throwNullArgument();
    }
}
