package pl.edu.mimuw.at429630.robson.conversion;

import lombok.Getter;
import lombok.Setter;
import lombok.ToString;
import lombok.experimental.Accessors;

import java.util.Objects;

@ToString
@Accessors(fluent = true)
public class VariableBuilder {
    @Getter
    @Setter
    private String name;

    public String make() {
        return String.format("public static Double %s = 0.0;", name);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        VariableBuilder that = (VariableBuilder) o;
        return that.name.equals(name);
    }

    @Override
    public int hashCode() {
        return Objects.hash(name);
    }
}
