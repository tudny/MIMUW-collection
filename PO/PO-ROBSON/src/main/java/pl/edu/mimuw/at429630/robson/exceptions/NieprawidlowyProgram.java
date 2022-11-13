package pl.edu.mimuw.at429630.robson.exceptions;

import java.io.IOException;

public class NieprawidlowyProgram extends Exception {
    public NieprawidlowyProgram(String message) {
        super(message);
    }

    public NieprawidlowyProgram(String s, Exception e) {
        super(s, e);
    }
}
