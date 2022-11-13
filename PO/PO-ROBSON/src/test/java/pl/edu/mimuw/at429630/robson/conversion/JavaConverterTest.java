package pl.edu.mimuw.at429630.robson.conversion;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import pl.edu.mimuw.at429630.robson.Robson;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;

import static org.junit.jupiter.api.Assertions.*;

class JavaConverterTest {

    @Test
    void TestConversion() throws NieprawidlowyProgram, BladWykonania {
        assertThrows(NieprawidlowyProgram.class, () -> {
            Robson robson = new Robson();
            robson.fromJSON("tests/official/empty.json");
        });

//        Double d = robson.wykonaj();
//        System.out.println(d);
//
//        JavaConverter converter = robson.toJava("");
//        System.out.println(converter);
//
//        System.out.println("\tpublic static Boolean javaDoubleToBoolean(Double value) {\n" +
//                "        if (value.equals(0.0))\n" +
//                "            return false;\n" +
//                "        else if (value.equals(1.0))\n" +
//                "            return true;\n" +
//                "\n" +
//                "        throw new RuntimeException(\"Trying to convert \" + value + \" to Boolean and this is not bool bro.\");\n" +
//                "    }\n" +
//                "\n" +
//                "    public static Double javaBooleanToDouble(Boolean value) {\n" +
//                "        return value ? 1.0 : 0.0;\n" +
//                "    }\n" +
//                "    \n" +
//                "    public static Double javaBooleanToDouble(Double value) {\n" +
//                "        return value;\n" +
//                "    }");
//
//        converter.getFunctions().forEach(functionBuilder -> {
//            System.out.println(functionBuilder.make());
//            System.out.println();
//        });
//
//        converter.getVariables().forEach(s -> {
//            System.out.println(s.make());
//        });
    }
}
