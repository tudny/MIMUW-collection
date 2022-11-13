package pl.edu.mimuw.at429630.robson.core;

import pl.edu.mimuw.at429630.robson.instructions.numbers.Number;

import java.util.HashMap;
import java.util.Map;

public class RuntimeContext {
    private final Map<String, Double> variables;

    public RuntimeContext() {
        variables = new HashMap<>();
    }

    public Double getVariable(String name) {
        return variables.getOrDefault(name, Number.DEFAULT_GLOBAL_VALUE);
    }

    public void setVariable(String name, Double value) {
        variables.put(name, value);
    }
}
