package concurrentcube.sequantial;

public class SideInfo {
    private final int opposite;
    private final int[] near;
    private final int[] howManyRotations;

    public SideInfo(int opposite, int[] near, int [] howManyRotations) {
        this.opposite = opposite;
        this.near = near;
        this.howManyRotations = howManyRotations;
    }

    public int getOpposite() {
        return opposite;
    }

    public int[] getNear() {
        return near;
    }

    public int[] getHowManyRotations() {
        return howManyRotations;
    }
}
