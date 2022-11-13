package pl.edu.mimuw.at429630.henry.robson;

public class LeftInstruction extends RobInstruction {
    @Override
    protected InstructionHandler getHandler() {
        return () -> controllerInstance.handleLeft();
    }
}
