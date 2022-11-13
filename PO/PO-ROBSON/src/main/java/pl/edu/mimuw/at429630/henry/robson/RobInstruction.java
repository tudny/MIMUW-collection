package pl.edu.mimuw.at429630.henry.robson;

import lombok.ToString;
import pl.edu.mimuw.at429630.henry.BoardController;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

@ToString
public abstract class RobInstruction extends Instruction {

    protected static BoardController controllerInstance;

    public static void setControllerInstance(BoardController controller) {
        controllerInstance = controller;
    }

    protected abstract InstructionHandler getHandler();

    @Override
    public Double execute(RuntimeContext context) throws BladWykonania {
        getHandler().handle();
        return 0.0;
    }

    @Override
    public String toJava(JavaConverter converter) {
        return "RobInstruction - not compilable.";
    }

    @Override
    public void validate() throws NieprawidlowyProgram { }
}
