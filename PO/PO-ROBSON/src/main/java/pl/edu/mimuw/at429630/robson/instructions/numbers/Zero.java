package pl.edu.mimuw.at429630.robson.instructions.numbers;

import lombok.ToString;

@ToString(callSuper = true)
public class Zero extends Number {
    public Zero() {
        super(0.0);
    }
}
