package pl.edu.mimuw.at429630.robson.instructions.numbers.comparators;

import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;

@FunctionalInterface
public interface NumberCompareOperation {
    Boolean operation(Double a, Double b) throws BladWykonania;
}
