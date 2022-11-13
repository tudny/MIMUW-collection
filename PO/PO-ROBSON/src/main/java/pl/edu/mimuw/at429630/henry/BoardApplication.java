package pl.edu.mimuw.at429630.henry;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Scene;
import javafx.scene.layout.Pane;
import javafx.stage.Stage;
import pl.edu.mimuw.at429630.henry.robson.RobInstruction;

public class BoardApplication extends Application {

    public static final String TITLE = "Henry the ROB";

    public static void run() {
        launch();
    }

    @Override
    public void start(Stage primaryStage) throws Exception {

        FXMLLoader loader = new FXMLLoader();
        loader.setLocation(getClass().getResource("/fxml/BoardWindow.fxml"));
        Pane root = loader.load();

        primaryStage.setTitle(TITLE);
        primaryStage.setScene(new Scene(root));
        primaryStage.setResizable(false);
        primaryStage.sizeToScene();
        primaryStage.setX(10.0);
        primaryStage.setY(10.0);

        BoardController controller = loader.getController();
        controller.setup(null, null, primaryStage);
        RobInstruction.setControllerInstance(controller);

        primaryStage.show();
    }
}
