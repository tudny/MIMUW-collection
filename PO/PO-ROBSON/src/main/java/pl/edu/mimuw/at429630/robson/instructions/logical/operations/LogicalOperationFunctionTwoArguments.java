package pl.edu.mimuw.at429630.robson.instructions.logical.operations;

import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;

@FunctionalInterface
public interface LogicalOperationFunctionTwoArguments {
    Double getOperation(Double a, Double b) throws BladWykonania;
}
