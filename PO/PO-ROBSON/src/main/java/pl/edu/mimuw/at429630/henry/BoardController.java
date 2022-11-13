package pl.edu.mimuw.at429630.henry;

import javafx.animation.AnimationTimer;
import javafx.collections.FXCollections;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.Node;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.Button;
import javafx.scene.control.TextField;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyEvent;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.RowConstraints;
import javafx.scene.text.Text;
import javafx.stage.Stage;
import pl.edu.mimuw.at429630.henry.robson.HenryRobsonConfig;
import pl.edu.mimuw.at429630.henry.robson.RobInstruction;
import pl.edu.mimuw.at429630.robson.Robson;
import pl.edu.mimuw.at429630.robson.exceptions.BladWykonania;
import pl.edu.mimuw.at429630.robson.exceptions.NieprawidlowyProgram;

import java.util.ArrayList;
import java.util.List;

public class BoardController {

    private static final Integer FOOD_ENERGY = 10;
    private static final Integer DEFAULT_X = 0;
    private static final Integer DEFAULT_Y = 0;
    private static final Integer INSTRUCTION_COST = 1;

    @FXML
    private GridPane boardGrid;

    @FXML
    private Text infoLabel;

    @FXML
    private GridPane footer;

    @FXML
    private Text energyInfo;

    /* COMMAND LINE BEGIN */

    @FXML
    private TextField commandLine;

    @FXML @SuppressWarnings("unused")
    private Button commitButton;

    @FXML
    void commitButtonPressed(@SuppressWarnings("unused") ActionEvent event) {
        newCommandTyped();
    }

    @FXML
    void keyPressedInCommandLine(KeyEvent event) {
        if (event.getCode() == KeyCode.ENTER) {
            newCommandTyped();
        }
    }

    private void newCommandTyped() {
        setInfoLabel(Info.INFO, String.format("Entered %s", commandLine.getText()));

        handleCommand(commandLine.getText());
        commandLine.clear();
    }

    /*
     * LOAD filepath
     * RUN
     * MOVE X Y
     * */
    private void handleCommand(String command) {

        try {
            String[] args = command.trim().split("\\s+");
            if (args.length < 1) {
                setInfoLabel(Info.INFO, "No command detected,");
            }

            switch (args[0]) {
                case "LOAD":
                    if (args.length != 2) {
                        setInfoLabel(Info.WARNING, "No file provided.");
                        break;
                    }
                    String filePath = args[1];
                    robson.fromJSON(filePath);
                    break;
                case "RUN":
                    Double d = robson.wykonaj();
                    System.out.println(d);
                    break;
                case "MOVE":
                    int x = Integer.parseInt(args[1]);
                    int y = Integer.parseInt(args[2]);
                    moveHenryTo(x, y);
                    break;
                default:
                    setInfoLabel(Info.WARNING, "Command " + args[0] + " is not known.");
            }

        } catch (IndexOutOfBoundsException | NieprawidlowyProgram | BladWykonania | NumberFormatException e) {
            setInfoLabel(Info.ERROR, e.getMessage());
        }
    }

    /* COMMAND LINE END */

    private interface HenryStep {
        void handle();
    }

    private final List<AnchorPane> fields = new ArrayList<>();
    private Integer width;
    private Integer height;

    private Stage primaryStage;

    private HenryView henry;
    private Image foodImage;

    private final List<ImageView> food = new ArrayList<>();

    private Robson robson;

    private AnimationTimer timer;
    private final List<HenryStep> steps = FXCollections.observableArrayList();

    private Double getWindowBorderWidth() {
        return primaryStage.getWidth() - primaryStage.getScene().getWidth();
    }

    private Double getWindowBorderHeight() {
        return primaryStage.getHeight() - primaryStage.getScene().getHeight() - footer.getHeight();
    }

    public enum Info {
        INFO,
        WARNING,
        ERROR
    }

    private void setInfoLabel(Info info, String message) {
        infoLabel.setText(String.format("[%s]: %s", info, message));
    }

    private void zeroAnchors(Node node) {
        AnchorPane.setTopAnchor(node, 0.0);
        AnchorPane.setBottomAnchor(node, 0.0);
        AnchorPane.setLeftAnchor(node, 0.0);
        AnchorPane.setRightAnchor(node, 0.0);
    }

    private void setupGrid(Integer inputWidth, Integer inputHeight) {

        if (inputWidth == null || inputHeight == null) {
            inputWidth = BoardConfig.DEFAULT_BOARD_WIDTH;
            inputHeight = BoardConfig.DEFAULT_BOARD_HEIGHT;
        }

        primaryStage.setWidth(inputWidth * BoardConfig.TAIL_WIDTH + getWindowBorderWidth());
        primaryStage.setHeight(inputHeight * BoardConfig.TAIL_HEIGHT + getWindowBorderHeight());

        double widthPercent = 100.0 / (double) inputWidth;
        double heightPercent = 100.0 / (double) inputHeight;

        boardGrid.getColumnConstraints().clear();
        boardGrid.getRowConstraints().clear();

        for (int i = 0; i < inputWidth; ++i) {
            ColumnConstraints columnConstraints = new ColumnConstraints();
            columnConstraints.setFillWidth(true);
            columnConstraints.setPercentWidth(widthPercent);
            boardGrid.getColumnConstraints().add(columnConstraints);
        }

        for (int i = 0; i < inputHeight; ++i) {
            RowConstraints rowConstraints = new RowConstraints();
            rowConstraints.setFillHeight(true);
            rowConstraints.setPercentHeight(heightPercent);
            boardGrid.getRowConstraints().add(rowConstraints);
        }

        Image terrain = new Image(getClass().getResourceAsStream("/textures/terrain_texture_border.jpg"));

        for (int x = 0; x < inputWidth; x++) {
            for (int y = 0; y < inputHeight; y++) {
                AnchorPane pane = new AnchorPane();
                ImageView tail = new ImageView(terrain);
                tail.setFitWidth(BoardConfig.TAIL_WIDTH);
                tail.setFitHeight(BoardConfig.TAIL_HEIGHT);
                pane.getChildren().add(tail);

                zeroAnchors(tail);

                boardGrid.add(pane, x, y);
                fields.add(pane);

                ImageView foodView = new ImageView(foodImage);
                food.add(foodView);
                pane.getChildren().add(foodView);
                zeroAnchors(foodView);
                foodView.setVisible(Boolean.FALSE);
                foodView.setFitWidth(BoardConfig.TAIL_WIDTH);
                foodView.setFitHeight(BoardConfig.TAIL_HEIGHT);

                final int _x = x, _y = y;
                tail.addEventHandler(MouseEvent.MOUSE_CLICKED, mouseEvent -> addFood(_x, _y));
            }
        }

        this.width = inputWidth;
        this.height = inputHeight;

        setInfoLabel(Info.INFO, String.format("Width: %d | Height: %d", boardGrid.getColumnCount(), boardGrid.getRowCount()));
    }

    private void setupHenry(int x, int y) {
        henry = new HenryView(x, y, 0, energyInfo, this::handleDead);
        moveHenryTo(x, y);
        RobInstruction.setControllerInstance(this);
    }

    private int getId(int x, int y) {
        return x * height + y;
    }

    private void moveHenryTo(int x, int y) {
        ImageView henryImageView = henry.getImageView();
        int lastId = getId(henry.getX(), henry.getY());
        int id = getId(x, y);

        if (!fields.get(id).getChildren().contains(henryImageView)) {
            fields.get(lastId).getChildren().removeAll(henryImageView);
            fields.get(id).getChildren().add(henryImageView);
        }

        if (hasFood(x, y)) {
            eatFood(x, y);
        }

        henry.setX(x);
        henry.setY(y);
    }

    private void setupFood() {
        foodImage = new Image(getClass().getResourceAsStream("/textures/battery.png"));
    }

    private void eatFood(int x, int y) {
        removeFood(x, y);
        henry.addEnergy(FOOD_ENERGY);
    }

    private Boolean hasFood(int x, int y) {
        return food.get(getId(x, y)).isVisible();
    }

    private void addFood(int x, int y) {
        food.get(getId(x, y)).setVisible(Boolean.TRUE);
        setInfoLabel(Info.INFO, "Added food on " + x + ", " + y + ".");
    }

    private void removeFood(int x, int y) {
        food.get(getId(x, y)).setVisible(Boolean.FALSE);
        setInfoLabel(Info.INFO, "Removed food on " + x + ", " + y + ".");
    }

    public void setup(Integer inputWidth, Integer inputHeight, Stage primaryStage) {
        this.primaryStage = primaryStage;
        this.robson = new Robson(new HenryRobsonConfig());

        setInfoLabel(Info.INFO, "Settings things up.");
        setupFood();
        setupGrid(inputWidth, inputHeight);
        setupHenry(DEFAULT_X, DEFAULT_Y);
        setupTimer();
    }

    private void setupTimer() {
        timer = new AnimationTimer() {
            public final long DELAY = 1_000_000_000L;
            private long lastFrame = -1;

            @Override
            public void handle(long l) {
                if (lastFrame == -1)
                    lastFrame = l;

                if (lastFrame + DELAY <= l) {
                    lastFrame = l;

                    if (!steps.isEmpty()) {
                        steps.remove(0).handle();
                        henry.takeEnergy(INSTRUCTION_COST);
                    }
                }
            }
        };

        timer.start();
    }

    public void handleLeft() {
        steps.add(() -> {
           System.out.println("left");
            setInfoLabel(Info.INFO, "Rotated to the left.");

            henry.rotate(-2);
        });
    }

    public void handleRight() {
        steps.add(() -> {
            System.out.println("right");
            setInfoLabel(Info.INFO, "Rotated to the right.");

            henry.rotate(2);
        });
    }

    public void handleForward() {
        steps.add(() -> {
            System.out.println("forward");
            setInfoLabel(Info.INFO, "Moved forward.");

            int x = (henry.getX() + henry.getDX() + width) % width;
            int y = (henry.getY() + henry.getDY() + height) % height;

            moveHenryTo(x, y);
        });
    }

    public void handleSmell() {
        steps.add(() -> {
            setInfoLabel(Info.INFO, "Smelled around");

            System.out.println("smell");
            for (int d = 0; d < 8; d += 2) {
                int x = (henry.getX() + HenryView.DX[d] + width) % width;
                int y = (henry.getY() + HenryView.DY[d] + height) % height;
                if (hasFood(x, y)) {
                    henry.setDirection(d);
                    break;
                }
            }
        });
    }

    public void handleEat() {
        steps.add(() -> {
            System.out.println("eat");
            setInfoLabel(Info.INFO, "Eaten the food.");

            for (int d = 0; d < HenryView.HENRY_IMAGES_COUNT; d++) {

                int x = (henry.getX() + HenryView.DX[d] + width) % width;
                int y = (henry.getY() + HenryView.DY[d] + height) % height;

                if (hasFood(x, y)) {
                    moveHenryTo(x, y);
                    break;
                }
            }
        });
    }

    private void handleDead() {
        timer.stop();
        timer = null;

        Alert alert = new Alert(AlertType.INFORMATION);
        alert.setTitle("Henry is over.");
        alert.setHeaderText("Henry is dead.");
        alert.setContentText("Unfortunately your instructions turned poor Henry into pice od dead DNA.");

        alert.setOnHidden(dialogEvent -> primaryStage.close());

        alert.show();
    }
}
