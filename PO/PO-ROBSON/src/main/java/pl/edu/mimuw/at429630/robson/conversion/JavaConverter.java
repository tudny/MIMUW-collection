package pl.edu.mimuw.at429630.robson.conversion;

import lombok.ToString;

import java.util.*;

@ToString
public class JavaConverter {

    private final Set<FunctionBuilder> functions;
    private final Set<VariableBuilder> variables;
    private Integer counter = 0;
    private static final String NAME = "Foo";

    public JavaConverter() {
        functions = new HashSet<>();
        variables = new HashSet<>();
    }

    public void addFunction(FunctionBuilder functionBuilder) {
        functions.add(functionBuilder);
    }

    public void addVariable(VariableBuilder variable) {
        variables.add(variable);
    }

    public Set<FunctionBuilder> getFunctions() {
        return new HashSet<>(functions);
    }

    public Set<VariableBuilder> getVariables() {
        return new HashSet<>(variables);
    }

    public String getNextName(Object obj) {
        return obj.getClass().getSimpleName() + NAME + counter++;
    }
}
