package zad1.at429630.ewolucja.config;

import zad1.at429630.ewolucja.board.Board;
import zad1.at429630.ewolucja.instr.Instruction;
import zad1.at429630.ewolucja.instr.Program;
import zad1.at429630.ewolucja.utils.Randomizer;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class ConfigBuilder {
    private final Config config = new Config();
    private Board board;

    public void readConfigFile(String pathToFile) throws FileNotFoundException {
        List<String> fileAsLines = readFile(pathToFile);

        for (String line : fileAsLines) {
            parseInput(line);
        }

        config.setRandomizer(new Randomizer());
    }

    public void readBoardFile(String pathToFile) throws FileNotFoundException {
        List<String> fileAsLines = readFile(pathToFile);

        int lineLength = fileAsLines.get(0).length();
        int lineCount = fileAsLines.size();

        board = new Board(lineLength, lineCount, fileAsLines);
        config.setRozmiarPlanszyX(board.getWidth());
        config.setRozmiarPlanszyY(board.getHeight());
    }

    private List<String> readFile(String pathToFile) throws FileNotFoundException {
        File file = new File(pathToFile);

        if (!file.exists())
            throw new FileNotFoundException("No such file " + pathToFile);

        List<String> fileAsLines = new ArrayList<>();
        try (Scanner scanner = new Scanner(file)) {
            while (scanner.hasNextLine()) {
                fileAsLines.add(scanner.nextLine());
            }
        }

        return fileAsLines;
    }

    private void parseInput(String line) {
        String[] splitLine = line.split(" ");
        String param = splitLine[0], argument;
        try {
            argument = splitLine[1];
        } catch (IndexOutOfBoundsException e) {
            argument = "";
        }

        switch (param) {
            case "ile_tur":
                config.setIleTur(Integer.valueOf(argument));
                break;
            case "pocz_ile_robów":
                config.setPoczIleRobow(Integer.valueOf(argument));
                break;
            case "pocz_progr":
                config.setPoczProg(Program.of(argument));
                break;
            case "pocz_energia":
                config.setPoczEnergia(Integer.valueOf(argument));
                break;
            case "ile_daje_jedzenie":
                config.setIleDajeJedzenie(Integer.valueOf(argument));
                break;
            case "ile_rośnie_jedzenie":
                config.setIleRosnieJedzenie(Integer.valueOf(argument));
                break;
            case "koszt_tury":
                config.setKosztTury(Integer.valueOf(argument));
                break;
            case "pr_powielenia":
                config.setPrPowielania(Double.valueOf(argument));
                break;
            case "ułamek_energii_rodzica":
                config.setUlamekEnergiiRodzica(Double.valueOf(argument));
                break;
            case "limit_powielania":
                config.setLimitPowielania(Integer.valueOf(argument));
                break;
            case "pr_usunięcia_instr":
                config.setPrUsunieciaInstr(Double.valueOf(argument));
                break;
            case "pr_dodania_instr":
                config.setPrDodanieInstr(Double.valueOf(argument));
                break;
            case "pr_zmiany_instr":
                config.setPrZmianaInstr(Double.valueOf(argument));
                break;
            case "spis_instr":
                config.setSpisInstr(Instruction.arrayOf(argument));
                break;
            case "co_ile_wypisz":
                config.setCoIleWypisz(Integer.valueOf(argument));
                break;
            default:
                throw new IllegalArgumentException("Passed argument is not a proper parameter!");
        }
    }

    public Config makeConfig() {
        return config;
    }

    public Board getBoard() {
        return board;
    }
}
