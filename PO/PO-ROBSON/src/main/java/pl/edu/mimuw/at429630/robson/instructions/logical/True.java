package pl.edu.mimuw.at429630.robson.instructions.logical;

import lombok.ToString;

@ToString(callSuper = true)
public class True extends Bool {
    True() {
        super(1.0);
    }
}
