package pl.edu.mimuw.at429630.henry.robson;

public class SmellInstruction extends RobInstruction {
    @Override
    protected InstructionHandler getHandler() {
        return () -> controllerInstance.handleSmell();
    }
}
