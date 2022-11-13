package pl.edu.mimuw.at429630.robson.instructions;

import lombok.ToString;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;

@ToString
public abstract class Instruction {
    public abstract Double execute(RuntimeContext context) throws BladWykonania;
    public abstract String toJava(JavaConverter converter);
    public abstract void validate() throws NieprawidlowyProgram;
    protected void throwNullArgument() throws NieprawidlowyProgram {
        throw new NieprawidlowyProgram("Null argument in " + getClass().getSimpleName());
    }
}
