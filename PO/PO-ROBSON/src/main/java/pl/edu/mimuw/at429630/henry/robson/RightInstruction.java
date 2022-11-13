package pl.edu.mimuw.at429630.henry.robson;

public class RightInstruction extends RobInstruction {
    @Override
    protected InstructionHandler getHandler() {
        return () -> controllerInstance.handleRight();
    }
}
