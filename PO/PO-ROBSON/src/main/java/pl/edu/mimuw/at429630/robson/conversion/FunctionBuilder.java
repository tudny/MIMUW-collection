package pl.edu.mimuw.at429630.robson.conversion;

import lombok.ToString;

import java.util.ArrayList;
import java.util.List;

@ToString
public class FunctionBuilder {

    private String name;
    private final List<String> lines = new ArrayList<>();
    private String toReturn;

    public void addLine(String line) {
        lines.add(line + ";");
    }

    public void addLineWithoutSemi(String line) {
        lines.add(line);
    }

    public void removeLastLineIfExists() {
        if (!lines.isEmpty())
            lines.remove(lines.size() - 1);
    }

    public void setToReturn(String line) {
        toReturn = line + ";";
    }

    public void setName(String name) {
        this.name = name;
    }

    public String make() {
        StringBuilder function = new StringBuilder();

        function.append(String.format("public static Double %s() {\n", name));
        lines.forEach(line -> {
            function.append(String.format("\t%s\n", line));
        });
        function.append(String.format("\treturn %s\n}\n", toReturn));

        return function.toString();
    }
}
