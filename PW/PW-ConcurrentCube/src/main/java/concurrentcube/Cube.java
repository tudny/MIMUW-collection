package concurrentcube;

import concurrentcube.sequantial.SequentialCube;

import java.util.concurrent.Semaphore;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.function.BiConsumer;

public class Cube {

    private final int size;
    private final BiConsumer<Integer, Integer> beforeRotation;
    private final BiConsumer<Integer, Integer> afterRotation;
    private final Runnable beforeShowing;
    private final Runnable afterShowing;
    private final SequentialCube cube;

    private final AtomicInteger currentlyRotatingGroup = new AtomicInteger(-1);
    private final AtomicInteger howManyRotates = new AtomicInteger(0);
    private final AtomicInteger howManyToExit = new AtomicInteger(0);

    private final Lock lock = new ReentrantLock(true);
    private final Condition entryCondition = lock.newCondition();
    private final Condition exitCondition = lock.newCondition();
    private final Semaphore[] layers;

    public Cube(int size,
                BiConsumer<Integer, Integer> beforeRotation,
                BiConsumer<Integer, Integer> afterRotation,
                Runnable beforeShowing, Runnable afterShowing) {
        this.size = size;
        this.beforeRotation = beforeRotation;
        this.afterRotation = afterRotation;
        this.beforeShowing = beforeShowing;
        this.afterShowing = afterShowing;

        this.cube = new SequentialCube(this.size);

        this.layers = new Semaphore[this.size];
        for (int i = 0; i < this.size; ++i) {
            this.layers[i] = new Semaphore(1);
        }
    }

    private int giveGroup(int side) {
        switch (side) {
            case 0:
            case 5:
                return 0;
            case 1:
            case 3:
                return 1;
            case 2:
            case 4:
                return 2;
        }
        return 3;
    }

    private int giveNormalizedLayer(int side, int layer) {
        return (side > SequentialCube.RELATION[side].getOpposite()) ? size - layer - 1 : layer;
    }

    private void entry(int id) throws InterruptedException {

        lock.lock();

        try {
            while (!(currentlyRotatingGroup.get() == id || currentlyRotatingGroup.get() == -1)) {
                entryCondition.await();
            }

            currentlyRotatingGroup.set(id);
            howManyRotates.incrementAndGet();
        }
        finally {
             lock.unlock();
        }
    }

    private void exit(int id) throws InterruptedException {
        lock.lock();

        try {
            howManyRotates.decrementAndGet();
            howManyToExit.incrementAndGet();

            while (howManyRotates.get() > 0) {
                exitCondition.await();
            }

            currentlyRotatingGroup.set(-2);
            exitCondition.signalAll();
        }
        finally {
            if (howManyToExit.decrementAndGet() == 0 && howManyRotates.get() == 0) {
                currentlyRotatingGroup.set(-1);
                entryCondition.signalAll();
            }

            lock.unlock();
        }

        if (Thread.currentThread().isInterrupted()) {
            throw new InterruptedException();
        }
    }

    public void rotate(int side, int layer) throws InterruptedException {
        int id = giveGroup(side);
        int forSem = giveNormalizedLayer(side, layer);

        entry(id);
        try {
            layers[forSem].acquire();
        } catch (InterruptedException e) {
            lock.lock();

            howManyRotates.getAndDecrement();
            if (howManyToExit.get() != 0 && howManyRotates.get() == 0) {
                currentlyRotatingGroup.set(-2);
                exitCondition.signalAll();
            }
            else if (howManyRotates.get() == 0 /* howManyToExit == 0 */) {
                currentlyRotatingGroup.set(-1);
                entryCondition.signalAll();
            }

            lock.unlock();
            throw e;
        }

        this.beforeRotation.accept(side, layer);
        this.cube.rotate(side, layer);
        this.afterRotation.accept(side, layer);

        layers[forSem].release();
        exit(id);
    }

    public String show() throws InterruptedException {
        int id = giveGroup(-1);

        entry(id);

        String toShow = "";
        try {
            this.beforeShowing.run();
            toShow = cube.show(); // add exception inside
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        } finally {
            this.afterShowing.run();
        }

        exit(id);

        return toShow;
    }
}
