package pl.edu.mimuw.at429630.henry.robson;

public class EatInstruction extends RobInstruction {
    @Override
    protected InstructionHandler getHandler() {
        return () -> controllerInstance.handleEat();
    }
}
