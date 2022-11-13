package pl.edu.mimuw.at429630.henry;

import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.AnchorPane;
import javafx.scene.text.Text;
import lombok.Getter;
import lombok.Setter;
import pl.edu.mimuw.at429630.henry.robson.DeadListener;

public class HenryView {

    public final static Integer[] DX = {0, 1, 1, 1, 0, -1, -1, -1};
    public final static Integer[] DY = {-1, -1, 0, 1, 1, 1, 0, -1};

    public final static Integer HENRY_IMAGES_COUNT = 8;
    private final Image[] henryImages = new Image[HENRY_IMAGES_COUNT];
    private final static Integer DEFAULT_ENERGY = 20;

    public Integer getDX() {
        return DX[direction];
    }

    public Integer getDY() {
        return DY[direction];
    }

    private void loadTextures() {
        for (int i = 0; i < HENRY_IMAGES_COUNT; i++) {
            String path = "/textures/henry_" + i + ".png";
            henryImages[i] = new Image(getClass().getResourceAsStream(path));
        }
    }

    @Getter
    private final ImageView imageView = new ImageView();

    @Getter
    private Integer direction;

    @Getter
    @Setter
    private Integer x, y;

    @Getter
    private int energy = 0;

    private final Text energyText;

    private final DeadListener deadListener;

    private void changeEnergy(int x) {
        energy += x;
        energyText.setText("Energy: " + energy);

        if (energy <= 0) {
            deadListener.handleDead();
        }
    }

    public void addEnergy(int x) {
        changeEnergy(x);
    }

    public void takeEnergy(int x) {
        changeEnergy(-x);
    }

    public HenryView(int x, int y, int direction, Text energyView, DeadListener deadListener) {
        loadTextures();

        imageView.setFitWidth(BoardConfig.TAIL_WIDTH);
        imageView.setFitHeight(BoardConfig.TAIL_HEIGHT);

        AnchorPane.setTopAnchor(imageView, 0.0);
        AnchorPane.setBottomAnchor(imageView, 0.0);
        AnchorPane.setLeftAnchor(imageView, 0.0);
        AnchorPane.setRightAnchor(imageView, 0.0);

        this.x = x;
        this.y = y;
        this.energyText = energyView;
        this.deadListener = deadListener;

        addEnergy(DEFAULT_ENERGY);

        setDirection(direction);
    }

    public void setDirection(Integer direction) {
        assert direction >= 0 && direction < HENRY_IMAGES_COUNT;

        this.direction = direction;

        updateDirectionImage();
    }

    public void rotate(int x) {
        setDirection((getDirection() + x + HENRY_IMAGES_COUNT) % HENRY_IMAGES_COUNT);
    }

    private void updateDirectionImage() {
        imageView.setImage(henryImages[direction]);
    }
}
