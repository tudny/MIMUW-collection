package pl.edu.mimuw.at429630.robson.core;

import com.squareup.moshi.adapters.PolymorphicJsonAdapterFactory;
import pl.edu.mimuw.at429630.robson.instructions.Block;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;
import pl.edu.mimuw.at429630.robson.instructions.logical.False;
import pl.edu.mimuw.at429630.robson.instructions.logical.If;
import pl.edu.mimuw.at429630.robson.instructions.logical.True;
import pl.edu.mimuw.at429630.robson.instructions.logical.operations.And;
import pl.edu.mimuw.at429630.robson.instructions.logical.operations.Not;
import pl.edu.mimuw.at429630.robson.instructions.logical.operations.Or;
import pl.edu.mimuw.at429630.robson.instructions.loops.While;
import pl.edu.mimuw.at429630.robson.instructions.numbers.Number;
import pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic.Divide;
import pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic.Minus;
import pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic.Plus;
import pl.edu.mimuw.at429630.robson.instructions.numbers.arithmetic.Times;
import pl.edu.mimuw.at429630.robson.instructions.numbers.comparators.*;
import pl.edu.mimuw.at429630.robson.instructions.variables.Assignment;
import pl.edu.mimuw.at429630.robson.instructions.variables.Variable;

public class InstructionJsonConfig {
    public InstructionJsonConfig() {}

    public PolymorphicJsonAdapterFactory<Instruction> getJsonAdapterFactory() {
        return PolymorphicJsonAdapterFactory.of(Instruction.class, "typ")
                /* Block */
                .withSubtype(Block.class, "Blok")
                /* Number (double) */
                .withSubtype(Number.class, "Liczba")
                /* Number operations */
                .withSubtype(Plus.class, "Plus")
                .withSubtype(Minus.class, "Minus")
                .withSubtype(Times.class, "Razy")
                .withSubtype(Divide.class, "Dzielenie")
                /* Booleans */
                .withSubtype(True.class, "True")
                .withSubtype(False.class, "False")
                /* Boolean operations */
                .withSubtype(And.class, "And")
                .withSubtype(Or.class, "Or")
                .withSubtype(Not.class,"Not")
                /* Number comparators */
                .withSubtype(LessThan.class, "<")
                .withSubtype(GreaterThan.class, ">")
                .withSubtype(LessOrEqual.class, "<=")
                .withSubtype(GreaterOrEqual.class, ">=")
                .withSubtype(Equals.class, "==")
                /* Loops */
                .withSubtype(While.class, "While")
                /* Variables */
                .withSubtype(Variable.class, "Zmienna")
                .withSubtype(Assignment.class, "Przypisanie")
                /* Logic */
                .withSubtype(If.class, "If");
    }
}
