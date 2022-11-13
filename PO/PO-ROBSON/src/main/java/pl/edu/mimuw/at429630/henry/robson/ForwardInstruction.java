package pl.edu.mimuw.at429630.henry.robson;

public class ForwardInstruction extends RobInstruction {
    @Override
    protected InstructionHandler getHandler() {
        return () -> controllerInstance.handleForward();
    }
}
