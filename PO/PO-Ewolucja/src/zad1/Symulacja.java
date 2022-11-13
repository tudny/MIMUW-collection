package zad1;

import zad1.at429630.ewolucja.Evolution;
import zad1.at429630.ewolucja.board.Board;
import zad1.at429630.ewolucja.config.Config;
import zad1.at429630.ewolucja.config.ConfigBuilder;

import java.io.FileNotFoundException;

/*
* Program wczytuje w plików parametry oraz planszę.
*
* zad1.Symulacja prowadzona jest w obiekcie klasy Evolution, gdzie znajduje się niemal cała logika.
* Roby mają dostęp do planszy, ponieważ poruszają się po niej (muszę widzieć co się dzieje wokół).
* W przypadku wkroczenia na nowe pole robiony jest callback, który informuje pole o wejściu Roba.
* */

public class Symulacja {
    public static void main(String[] args) {

        if (args.length != 2) {
            System.err.println("The execution should be: java zad1/zad1.Symulacja plansza.txt parametry.txt");
            System.err.printf("Provided %d arguments instead of 2.", args.length);
            System.exit(1);
        }

        try {
            String boardFilePath = args[0];
            String parametersFilePath = args[1];

            ConfigBuilder configBuilder = new ConfigBuilder();
            configBuilder.readConfigFile(parametersFilePath);
            configBuilder.readBoardFile(boardFilePath);

            Config evolutionConfig = configBuilder.makeConfig();
            Board board = configBuilder.getBoard();

            Evolution evolution = new Evolution(evolutionConfig, board);
            evolution.fillRobs();
            evolution.evolve();
        } catch (FileNotFoundException exception) {
            System.err.println("Provided files don't exist. " + exception.getMessage());
            System.exit(2);
        }
    }
}
