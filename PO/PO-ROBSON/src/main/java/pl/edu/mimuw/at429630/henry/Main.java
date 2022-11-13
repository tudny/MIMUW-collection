package pl.edu.mimuw.at429630.henry;

public class Main {
    public static void main(String[] args) {
        for (String arg : args) {
            try {
                int value = Integer.parseInt(arg.substring(3));

                if (arg.startsWith("-TW")) {
                    BoardConfig.TAIL_WIDTH = value;
                }
                else if (arg.startsWith("-TH")) {
                    BoardConfig.TAIL_HEIGHT = value;
                }
                else if (arg.startsWith("-BW")) {
                    BoardConfig.DEFAULT_BOARD_WIDTH = value;
                }
                else if (arg.startsWith("-BH")) {
                    BoardConfig.DEFAULT_BOARD_HEIGHT = value;
                }
            } catch (NumberFormatException e) {
                System.err.println("Program runs all good, but");
                System.err.println(e.getMessage());
            }
        }

        BoardApplication.run();
    }
}
