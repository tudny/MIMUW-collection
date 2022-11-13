package pl.edu.mimuw.at429630.robson.instructions.logical.operations;

import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;

public interface LogicalOperationFunctionSingleArgument {
    Double getOperation(Double a) throws BladWykonania;
}
