package pl.edu.mimuw.at429630.robson;

import com.squareup.moshi.JsonAdapter;
import com.squareup.moshi.JsonDataException;
import com.squareup.moshi.Moshi;
import pl.edu.mimuw.at429630.robson.conversion.JavaConverter;
import pl.edu.mimuw.at429630.robson.core.InstructionJsonConfig;
import pl.edu.mimuw.at429630.robson.core.RuntimeContext;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;
import pl.edu.mimuw.at429630.robson.instructions.Instruction;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Objects;
import java.util.Scanner;
import java.util.stream.Collectors;

public class Robson {
    private Instruction program;
    private final JsonAdapter<Instruction> jsonAdapter;

    public Robson() {
        this(new InstructionJsonConfig());
    }

    public Robson(InstructionJsonConfig config) {
        Moshi moshi = new Moshi.Builder()
                .add(config.getJsonAdapterFactory())
                .build();

        jsonAdapter = moshi.adapter(Instruction.class);
    }

    public void fromJSON(String filename)throws NieprawidlowyProgram {
        File file = new File(filename);
        try (InputStream in = new FileInputStream(file)) {
            String json = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))
                    .lines()
                    .collect(Collectors.joining("\n"));

            try {
                program = jsonAdapter.fromJson(json);
                Objects.requireNonNull(program).validate();
            } catch (IOException | JsonDataException | NullPointerException | NieprawidlowyProgram e) {
                throw new NieprawidlowyProgram("Error while parsing program.", e);
            }
        } catch (IOException e) {
            throw new NieprawidlowyProgram("File is not found.", e);
        }
    }

    public void toJSON(String filename) {
        File file = new File(filename);
        try {
            if (!file.createNewFile()) {
                System.out.println("Overriding file.");
            }

            String json = jsonAdapter.toJson(program);
            try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(file))) {
                bufferedWriter.write(json);
            }
        } catch (IOException e) {
            System.err.println(e.getMessage() + " with " + file.getAbsolutePath());
            e.printStackTrace();
        }
    }

    public double wykonaj() throws BladWykonania {
        if (program == null)
            throw new BladWykonania("No program was loaded.", new NullPointerException("Null program."));

        RuntimeContext context = new RuntimeContext();
        return program.execute(context);
    }

    public void toJava(String filename) {
        JavaConverter converter = new JavaConverter();
        String startPoint = program.toJava(converter);

        try (InputStream inputStream = getClass().getResourceAsStream("/template.java")) {
            String template = new BufferedReader(new InputStreamReader(inputStream, StandardCharsets.UTF_8))
                    .lines()
                    .collect(Collectors.joining("\n"))
                    .replace("%function_name%",startPoint + "()");

            File javaFile = new File(filename);
            if (!javaFile.createNewFile()) {
                System.out.println("File will be overridden.");
            }

            StringBuilder classBodyBuilder = new StringBuilder(template + "\n\n");

            converter.getVariables().forEach(variableBuilder -> {
                classBodyBuilder.append(variableBuilder.make()).append("\n");
            });

            converter.getFunctions().forEach(functionBuilder -> {
                classBodyBuilder.append(functionBuilder.make()).append("\n");
            });

            String classBody = classBodyBuilder.toString();

            StringBuilder wholeClassBuilder = new StringBuilder();

            wholeClassBuilder.append(String.format("class %s {\n", filename.replace(".java", "")));

            Scanner reader = new Scanner(classBody);
            while (reader.hasNextLine()) {
                String line = "\t" + reader.nextLine() + "\n";
                wholeClassBuilder.append(line);
            }

            wholeClassBuilder.append("}");

            String wholeClass = wholeClassBuilder.toString();

            try (PrintWriter printWriter = new PrintWriter(javaFile)) {
                printWriter.write(wholeClass);
                printWriter.flush();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
