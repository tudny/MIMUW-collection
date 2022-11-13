package concurrentcube;

import concurrentcube.sequantial.SequentialCube;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.function.Executable;

import java.util.ArrayList;
import java.util.Objects;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.junit.jupiter.api.Assertions.*;

public class CubeTest {

    private static String numberToColor(char c) {
        String ret = "";
        switch (c) {
            case '0':
                ret += "\u001b[37m";
                break;
            case '1':
                ret += "\u001b[31m";
                break;
            case '2':
                ret += "\u001b[34m";
                break;
            case '3':
                ret += "\u001b[35m";
                break;
            case '4':
                ret += "\u001b[32m";
                break;
            case '5':
                ret += "\u001b[33m";
                break;
        }

        return ret + c + "\u001b[0m";
    }

    public static boolean countCubeColor(String cube, int n) {
        return cube.chars().mapToObj(e -> (char) e)
                .collect(Collectors.groupingBy(Function.identity(), Collectors.counting()))
                .entrySet().stream()
                .filter(charCount -> charCount.getValue() == (long) n * n)
                .count() == 6;
    }

    public static Cube generateSimpleCube(int size) {
        return new Cube(size, (a, b) -> {}, (a, b) -> {}, () -> {}, () -> {});
    }

    public static String showColored(String cube) {
        int n = (int) Math.sqrt((double)(cube.length() / 6));

        StringBuilder colored = new StringBuilder();
        for (int i = 0; i < cube.length(); ++i) {
            if (i % n == 0) colored.append("\n");
            if (i % (n * n) == 0) colored.append("\n");
            colored.append(numberToColor(cube.charAt(i)));
        }

        return colored.toString();
    }

    @Test
    public void testThreads() {
        Cube cube = new Cube(3,
                (side, layer) -> { },
                (side, layer) -> { },
                () -> { },
                () -> { }
        );

        Thread t1 = new Thread(() -> {
            try {
                cube.rotate(0, 0);
                cube.rotate(1, 0);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }, "t1");

        Thread t2 = new Thread(() -> {
            try {
                cube.rotate(0, 1);
                cube.rotate(1, 1);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }, "t2");

        Thread t3 = new Thread(() -> {
            try {
                cube.rotate(1, 0);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }, "t3");

        t3.start();
        t1.start();
        t2.start();

        try {
            t1.join();
            t2.join();
            t3.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        try {
            assertTrue(countCubeColor(cube.show(), 3));
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    @Test
    @DisplayName("Simple Sequential Test")
    public void testSimpleSequentialRotate() {
        int size = 10;
        Cube cube = generateSimpleCube(size);

        Thread t1 = new Thread(() -> {
            for (int side = 0; side < 6; ++side) {
                for (int layer = 0; layer < size; layer++) {
                    try {
                        cube.rotate(side, layer);
                    } catch (InterruptedException ignored) {}
                }
            }
        }, "t1");

        t1.start();
        assertDoesNotThrow((Executable) t1::join);
        final String[] cubeShow = { "" };
        assertDoesNotThrow(() -> { cubeShow[0] = cube.show(); });

        assertTrue(countCubeColor(cubeShow[0], size));
    }

    @Test
    @DisplayName("Simple Sequential Test 2")
    public void testSimpleSequentialRotate2() {
        int size = 3;
        Cube cube = generateSimpleCube(size);

        Thread t1 = new Thread(() -> {
            try {
                cube.rotate(0, 1);
                cube.rotate(1, 2);
                cube.rotate(2, 0);
                cube.rotate(3, 1);
                cube.rotate(4, 2);
                cube.rotate(5, 0);
            } catch (InterruptedException ignored) {}
        }, "t1");

        t1.start();
        assertDoesNotThrow((Executable) t1::join);
        final String[] cubeShow = { "" };
        assertDoesNotThrow(() -> { cubeShow[0] = cube.show(); });

        assertTrue(countCubeColor(cubeShow[0], size));
        assertEquals("034031004111220504220355111343443220524501343555415232", cubeShow[0]);
    }

    @Test
    @DisplayName("Simple Concurrent Test")
    public void testConcurrentTest() {
        final int size = 200;
        final int threadCount = 1000;
        AtomicInteger beforeCounter = new AtomicInteger();
        AtomicInteger afterCounter = new AtomicInteger();
        Cube cube = new Cube(size,
                (s, l) -> beforeCounter.incrementAndGet(),
                (s, l) -> afterCounter.incrementAndGet(),
                () -> {},
                () -> {}
        );

        Thread[] threads = new Thread[threadCount];
        for (int i = 0; i < threadCount; ++i) {
            int finalI = i;
            threads[i] = new Thread(() -> {
                try {
                    cube.rotate(0, finalI % size);
                } catch (InterruptedException ignored) {}
            }, "t" + i);
        }

        for (Thread thread : threads) {
            thread.start();
        }

        try {
            for (Thread thread : threads) {
                thread.join();
            }
        } catch (InterruptedException ignored) {}

        final String[] cubeShow = new String[1];
        assertDoesNotThrow(() -> { cubeShow[0] = cube.show(); });
        assertTrue(countCubeColor(cubeShow[0], size));

        assertEquals(threadCount, beforeCounter.get());
        assertEquals(threadCount, afterCounter.get());
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

    private int giveNormalizedLayer(int side, int layer, int size) {
        return (side > SequentialCube.RELATION[side].getOpposite()) ? size - layer - 1 : layer;
    }

    @Test
    @DisplayName("Simple Concurrent Test 2")
    public void testConcurrentTest2() {
        final int size = 3;
        final int threadCount = 1000;
        AtomicInteger beforeCounter = new AtomicInteger();
        AtomicInteger afterCounter = new AtomicInteger();
        final AtomicReferenceArray<Integer> inside = new AtomicReferenceArray<>(new Integer[size]);
        final ConcurrentLinkedQueue<Integer> insideRapport = new ConcurrentLinkedQueue<>();
        final ConcurrentHashMap<Integer, AtomicInteger> sider = new ConcurrentHashMap<>();
        final Boolean[] error = { false };
        final Semaphore sem = new Semaphore(1);

        for (int i = 0; i < size; i++) {
            inside.set(i, 0);
        }

        for (int i = 0; i < 6; i++) {
            sider.put(i, new AtomicInteger());
        }

        Cube cube = new Cube(size,
                (s, l) -> {
                            int id = giveGroup(s);
                            int forSem = giveNormalizedLayer(s, l, size);
                            beforeCounter.incrementAndGet();
                            int ins = inside.getAndUpdate(forSem, integer -> integer + 1);
                            insideRapport.add(ins);
                            try {
                                sem.acquire();
                                sider.get(id).incrementAndGet();
                                for (int i = 0; i < 6; i++) {
                                    if (i != id && sider.get(i).get() != 0) {
                                        error[0] = true;
                                    }
                                }
                            } catch (InterruptedException e) {
                                e.printStackTrace();
                            } finally {
                                sem.release();
                            }
                        },
                (s, l) -> {
                            int id = giveGroup(s);
                            int forSem = giveNormalizedLayer(s, l, size);
                            afterCounter.incrementAndGet();
                            int ins = inside.updateAndGet(forSem, integer -> integer - 1);
                            insideRapport.add(ins);
                            try {
                                sem.acquire();
                                sider.get(id).decrementAndGet();
                                for (int i = 0; i < 4; i++) {
                                    if (i != id && sider.get(i).get() != 0) {
                                        System.out.println("Im " + s + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                        error[0] = true;
                                    }
                                }
                            } catch (InterruptedException e) {
                                e.printStackTrace();
                            } finally {
                                sem.release();
                            }
                        },
                () -> {
                            int id = giveGroup(-1);
                            try {
                                sem.acquire();
                                sider.get(id).incrementAndGet();
                                for (int i = 0; i < 6; i++) {
                                    if (i != id && sider.get(i).get() != 0) {
                                        error[0] = true;
                                    }
                                }
                            } catch (InterruptedException e) {
                                e.printStackTrace();
                            } finally {
                                sem.release();
                            }
                },
                () -> {
                            int id = giveGroup(-1);
                            try {
                                sem.acquire();
                                sider.get(id).decrementAndGet();
                                for (int i = 0; i < 4; i++) {
                                    if (i != id && sider.get(i).get() != 0) {
                                        System.out.println("Im " + -1 + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                        error[0] = true;
                                    }
                                }
                            } catch (InterruptedException e) {
                                e.printStackTrace();
                            } finally {
                                sem.release();
                            }
                }
        );

        Thread[] threads = new Thread[threadCount];
        for (int i = 0; i < threadCount; ++i) {
            int finalI = i;
            threads[i] = new Thread(() -> {
                try {
                    cube.rotate(finalI % 6, finalI % size);
                    String unused = cube.show();
                } catch (InterruptedException ignored) {}
            }, "t" + i);
        }

        for (Thread thread : threads) {
            thread.start();
        }

        try {
            for (Thread thread : threads) {
                thread.join();
            }
        } catch (InterruptedException ignored) {}

        final String[] cubeShow = new String[1];
        assertDoesNotThrow(() -> { cubeShow[0] = cube.show(); });
        assertTrue(countCubeColor(cubeShow[0], size));

        assertEquals(threadCount, beforeCounter.get());
        assertEquals(threadCount, afterCounter.get());

        var rapport = insideRapport.stream().collect(Collectors.groupingBy( Function.identity(), Collectors.counting()));
        assertTrue(rapport.containsKey(0) && rapport.size() == 1);

        assertFalse(error[0]);
    }

    /* http://www.java2s.com/example/java-utility-method/integer-array-permutation/nextpermutation-int-next-14492.html */
    /* Funkcja analogiczna do funkcji std::next_permutation z C++ */
    public static int[] nextPermutation(int[] next) {
        // the counts can swap among each other. The int[] is originally in ascending order
        // this generates the next array in lexicographic order descending

        // locate the last occurrence where next[k] < next[k+1]
        int gt = -1;
        for (int idx = 0; idx < next.length - 1; idx++) {
            if (next[idx] < next[idx + 1]) {
                gt = idx;//from  w  w  w  . j  a v a  2 s  .c om
            }
        }

        if (gt == -1) {
            return null;
        }

        int largestLessThan = gt + 1;
        for (int idx = 1 + largestLessThan; idx < next.length; idx++) {
            if (next[gt] < next[idx]) {
                largestLessThan = idx;
            }
        }

        int val = next[gt];
        next[gt] = next[largestLessThan];
        next[largestLessThan] = val;

        // reverse the tail of the array
        int[] newTail = new int[next.length - gt - 1];
        int ctr = 0;
        for (int idx = next.length - 1; idx > gt; idx--) {
            newTail[ctr++] = next[idx];
        }

        for (int idx = 0; idx < newTail.length; idx++) {
            next[gt + idx + 1] = newTail[idx];
        }

        return next;
    }

    private static class Move {
        private final int side;
        private final int layer;

        public Move(int side, int layer) {
            this.side = side;
            this.layer = layer;
        }

        public int getSide() {
            return side;
        }

        public int getLayer() {
            return layer;
        }

        public void apply(Cube cube) throws InterruptedException {
            cube.rotate(side, layer);
        }

        @Override
        public String toString() {
            return "(" +
                    "side=" + side +
                    ", layer=" + layer +
                    ')';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Move move = (Move) o;
            return side == move.side && layer == move.layer;
        }

        @Override
        public int hashCode() {
            return Objects.hash(side, layer);
        }

        public boolean beforeOfEqualByLayer(Move move) {
            return this.layer <= move.layer;
        }
    }

    @Test
    @DisplayName("Concurrent Scenario")
    public void testConcurrentScenario() {
        int nrOfMoves = 4;
        int size = 3;
        final ArrayList<Move> moves = new ArrayList<>();
        for (int i = 0; i < nrOfMoves; ++i) {
            moves.add(new Move(i % 6, i % size));
            System.out.println(moves.get(i));
        }

        final Cube cube1 = generateSimpleCube(size);

        Thread[] rotators = new Thread[nrOfMoves];
        for (int i = 0; i < nrOfMoves; ++i) {
            int finalI = i;
            rotators[i] = new Thread(() -> {
                try {
                    moves.get(finalI).apply(cube1);
                } catch (InterruptedException ignored) {}
            });
        }

        for (Thread t : rotators) {
            t.start();
        }

        for (Thread t : rotators) {
            try {
                t.join();
            } catch (InterruptedException ignored) {}
        }

        String cube1Show = "";
        try {
             cube1Show = cube1.show();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        assertTrue(countCubeColor(cube1Show, size));

        System.out.println(cube1Show);

        int[] permutation = IntStream.range(0, nrOfMoves).toArray();

        for (int i = 0; i < factorial(nrOfMoves); i++) {
            Cube cube = generateSimpleCube(size);

            for (int m = 0; m < permutation.length; m++) {
                try {
                    moves.get(permutation[m]).apply(cube);
                } catch (InterruptedException ignored) {}
            }

            try {
                String shown = cube.show();
                System.out.println(shown);

                if (shown.equals(cube1Show)) {
                    return;
                }

            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            permutation = nextPermutation(permutation);
        }

        fail();
    }

    private long factorial(int nrOfMoves) {
        if (nrOfMoves < 2) return 1;
        return factorial(nrOfMoves - 1) * nrOfMoves;
    }

    @Test
    @DisplayName("Concurrent Test with Interruptions")
    public void testConcurrentInterruptionsTest() {
        final int size = 20;
        final int threadCount = 10000;
        AtomicInteger beforeCounter = new AtomicInteger();
        AtomicInteger afterCounter = new AtomicInteger();
        final AtomicReferenceArray<Integer> inside = new AtomicReferenceArray<>(new Integer[size]);
        final ConcurrentLinkedQueue<Integer> insideRapport = new ConcurrentLinkedQueue<>();
        final ConcurrentHashMap<Integer, AtomicInteger> sider = new ConcurrentHashMap<>();
        final Boolean[] error = { false };
        final Semaphore sem = new Semaphore(1);

        for (int i = 0; i < size; i++) {
            inside.set(i, 0);
        }

        for (int i = 0; i < 6; i++) {
            sider.put(i, new AtomicInteger());
        }

        Cube cube = new Cube(size,
                (s, l) -> {
                    int id = giveGroup(s);
                    int forSem = giveNormalizedLayer(s, l, size);
                    beforeCounter.incrementAndGet();
                    int ins = inside.getAndUpdate(forSem, integer -> integer + 1);
                    insideRapport.add(ins);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).incrementAndGet();
                        for (int i = 0; i < 6; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                (s, l) -> {
                    int id = giveGroup(s);
                    int forSem = giveNormalizedLayer(s, l, size);
                    afterCounter.incrementAndGet();
                    int ins = inside.updateAndGet(forSem, integer -> integer - 1);
                    insideRapport.add(ins);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).decrementAndGet();
                        for (int i = 0; i < 4; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                System.out.println("Im " + s + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                () -> {
                    int id = giveGroup(-1);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).incrementAndGet();
                        for (int i = 0; i < 6; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                () -> {
                    int id = giveGroup(-1);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).decrementAndGet();
                        for (int i = 0; i < 4; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                System.out.println("Im " + -1 + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                }
        );

        Thread[] threads = new Thread[threadCount];
        for (int i = 0; i < threadCount; ++i) {
            int finalI = i;
            threads[i] = new Thread(() -> {
                try {
                    cube.rotate(finalI % 6, finalI % size);
                    String unused = cube.show();
                } catch (InterruptedException ignored) {}
            }, "t" + i);
        }

        for (Thread thread : threads) {
            thread.start();
        }

        try {
            Thread.sleep(1000);
        } catch (InterruptedException ignored) {}

        for (int i = 0; i < threads.length; i += 2) {
            threads[i].interrupt();
        }

        try {
            for (Thread thread : threads) {
                thread.join();
            }
        } catch (InterruptedException ignored) {}

        final String[] cubeShow = new String[1];
        assertDoesNotThrow(() -> { cubeShow[0] = cube.show(); });
        assertTrue(countCubeColor(cubeShow[0], size));

        assertEquals(afterCounter.get(), beforeCounter.get());

        var rapport = insideRapport.stream().collect(Collectors.groupingBy( Function.identity(), Collectors.counting()));
        assertTrue(rapport.containsKey(0) && rapport.size() == 1);

        assertFalse(error[0]);
    }

    @Test
    @DisplayName("Concurrent Test with Interruptions 2")
    public void testConcurrentInterruptions2Test() {
        final int size = 20;
        final int threadCount = 10000;
        AtomicInteger beforeCounter = new AtomicInteger();
        AtomicInteger afterCounter = new AtomicInteger();
        final AtomicReferenceArray<Integer> inside = new AtomicReferenceArray<>(new Integer[size]);
        final ConcurrentLinkedQueue<Integer> insideRapport = new ConcurrentLinkedQueue<>();
        final ConcurrentHashMap<Integer, AtomicInteger> sider = new ConcurrentHashMap<>();
        final Boolean[] error = { false };
        final Semaphore sem = new Semaphore(1);

        for (int i = 0; i < size; i++) {
            inside.set(i, 0);
        }

        for (int i = 0; i < 6; i++) {
            sider.put(i, new AtomicInteger());
        }

        Cube cube = new Cube(size,
                (s, l) -> {
                    int id = giveGroup(s);
                    int forSem = giveNormalizedLayer(s, l, size);
                    beforeCounter.incrementAndGet();
                    int ins = inside.getAndUpdate(forSem, integer -> integer + 1);
                    insideRapport.add(ins);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).incrementAndGet();
                        for (int i = 0; i < 6; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                (s, l) -> {
                    int id = giveGroup(s);
                    int forSem = giveNormalizedLayer(s, l, size);
                    afterCounter.incrementAndGet();
                    int ins = inside.updateAndGet(forSem, integer -> integer - 1);
                    insideRapport.add(ins);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).decrementAndGet();
                        for (int i = 0; i < 4; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                System.out.println("Im " + s + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                () -> {
                    int id = giveGroup(-1);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).incrementAndGet();
                        for (int i = 0; i < 6; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                () -> {
                    int id = giveGroup(-1);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).decrementAndGet();
                        for (int i = 0; i < 4; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                System.out.println("Im " + -1 + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                }
        );

        Thread[] threads = new Thread[threadCount];
        for (int i = 0; i < threadCount; ++i) {
            int finalI = i;
            threads[i] = new Thread(() -> {
                try {
                    cube.rotate(finalI % 6, finalI % size);
                    String unused = cube.show();
                } catch (InterruptedException ignored) {}
            }, "t" + i);
        }

        for (Thread thread : threads) {
            thread.start();
        }

        try {
            Thread.sleep(1000);
        } catch (InterruptedException ignored) {}

        for (Thread value : threads) {
            value.interrupt();
        }

        try {
            for (Thread thread : threads) {
                thread.join();
            }
        } catch (InterruptedException ignored) {}

        final String[] cubeShow = new String[1];
        assertDoesNotThrow(() -> { cubeShow[0] = cube.show(); });
        assertTrue(countCubeColor(cubeShow[0], size));

        assertEquals(afterCounter.get(), beforeCounter.get());

        var rapport = insideRapport.stream().collect(Collectors.groupingBy( Function.identity(), Collectors.counting()));
        assertTrue(rapport.containsKey(0) && rapport.size() == 1);

        assertFalse(error[0]);
    }


    @Test
    @DisplayName("Multiple Rotations in Threads")
    public void testMultipleRotationsInThreads() {
        final int size = 10;
        final int threadCount = 1000;
        final int operationPerThread = 2000;

        AtomicInteger beforeCounter = new AtomicInteger();
        AtomicInteger afterCounter = new AtomicInteger();
        final AtomicReferenceArray<Integer> inside = new AtomicReferenceArray<>(new Integer[size]);
        final ConcurrentLinkedQueue<Integer> insideRapport = new ConcurrentLinkedQueue<>();
        final ConcurrentHashMap<Integer, AtomicInteger> sider = new ConcurrentHashMap<>();
        final Boolean[] error = { false };
        final Semaphore sem = new Semaphore(1);

        for (int i = 0; i < size; i++) {
            inside.set(i, 0);
        }

        for (int i = 0; i < 6; i++) {
            sider.put(i, new AtomicInteger());
        }

        Cube cube = new Cube(size,
                (s, l) -> {
                    int id = giveGroup(s);
                    int forSem = giveNormalizedLayer(s, l, size);
                    beforeCounter.incrementAndGet();
                    int ins = inside.getAndUpdate(forSem, integer -> integer + 1);
                    insideRapport.add(ins);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).incrementAndGet();
                        for (int i = 0; i < 6; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                (s, l) -> {
                    int id = giveGroup(s);
                    int forSem = giveNormalizedLayer(s, l, size);
                    afterCounter.incrementAndGet();
                    int ins = inside.updateAndGet(forSem, integer -> integer - 1);
                    insideRapport.add(ins);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).decrementAndGet();
                        for (int i = 0; i < 4; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                System.out.println("Im " + s + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                () -> {
                    int id = giveGroup(-1);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).incrementAndGet();
                        for (int i = 0; i < 6; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                () -> {
                    int id = giveGroup(-1);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).decrementAndGet();
                        for (int i = 0; i < 4; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                System.out.println("Im " + -1 + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                }
        );

        Thread[] threads = new Thread[threadCount];
        for (int i = 0; i < threadCount; ++i) {
            threads[i] = new Thread(() -> {
                for (int o = 0; o < operationPerThread; ++o) {
                    try {
                        cube.rotate((1324 * o) % 6, (738 * o) % size);
                    } catch (InterruptedException ignored) {}
                }
            }, "t" + i);
        }

        for (int i = 0; i < threadCount; ++i) {
            threads[i].start();
        }

        try {
            for (int i = 0; i < threadCount; ++i) {
                threads[i].join();
            }
        } catch (InterruptedException ignored) {}


        final String[] cubeShow = new String[1];
        assertDoesNotThrow(() -> { cubeShow[0] = cube.show(); });
        assertTrue(countCubeColor(cubeShow[0], size));

        assertEquals(afterCounter.get(), beforeCounter.get());

        var rapport = insideRapport.stream().collect(Collectors.groupingBy( Function.identity(), Collectors.counting()));
        assertTrue(rapport.containsKey(0) && rapport.size() == 1);

        assertFalse(error[0]);
    }

    @Test
    @DisplayName("Not On Blocking Barrier")
    public void testNotBlockingBarrier() {
        final int size = 1000;
        final int threadCount = 5;

        AtomicInteger beforeCounter = new AtomicInteger();
        AtomicInteger afterCounter = new AtomicInteger();
        final AtomicReferenceArray<Integer> inside = new AtomicReferenceArray<>(new Integer[size]);
        final ConcurrentLinkedQueue<Integer> insideRapport = new ConcurrentLinkedQueue<>();
        final ConcurrentHashMap<Integer, AtomicInteger> sider = new ConcurrentHashMap<>();
        final Boolean[] error = { false };
        final Semaphore sem = new Semaphore(1);
        final ConcurrentLinkedQueue<Move> logger = new ConcurrentLinkedQueue<>();

        for (int i = 0; i < size; i++) {
            inside.set(i, 0);
        }

        for (int i = 0; i < 6; i++) {
            sider.put(i, new AtomicInteger());
        }

        Cube cube = new Cube(size,
                (s, l) -> {
                    logger.add(new Move(s, l));
                    int id = giveGroup(s);
                    int forSem = giveNormalizedLayer(s, l, size);
                    beforeCounter.incrementAndGet();
                    int ins = inside.getAndUpdate(forSem, integer -> integer + 1);
                    insideRapport.add(ins);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).incrementAndGet();
                        for (int i = 0; i < 6; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                (s, l) -> {
                    logger.add(new Move(s, l));
                    int id = giveGroup(s);
                    int forSem = giveNormalizedLayer(s, l, size);
                    afterCounter.incrementAndGet();
                    int ins = inside.updateAndGet(forSem, integer -> integer - 1);
                    insideRapport.add(ins);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).decrementAndGet();
                        for (int i = 0; i < 4; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                System.out.println("Im " + s + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                () -> {
                    int id = giveGroup(-1);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).incrementAndGet();
                        for (int i = 0; i < 6; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                },
                () -> {
                    int id = giveGroup(-1);
                    try {
                        sem.acquireUninterruptibly();
                        sider.get(id).decrementAndGet();
                        for (int i = 0; i < 4; i++) {
                            if (i != id && sider.get(i).get() != 0) {
                                System.out.println("Im " + -1 + "(" + id + ") and " + i + " has " + sider.get(i).get());
                                error[0] = true;
                            }
                        }
                    } finally {
                        sem.release();
                    }
                }
        );

        Thread[] threads = new Thread[threadCount];
        Lock lock = new ReentrantLock();
        Condition condition = lock.newCondition();
        final AtomicInteger whoRotates = new AtomicInteger(threadCount - 1);

        for (int i = 0; i < threadCount; ++i) {
            int finalI = i;
            threads[i] = new Thread(() -> {
                try {
                    lock.lock();
                    while (whoRotates.get() != finalI) {
                        condition.await();
                    }
                    cube.rotate(0, finalI);

                    whoRotates.decrementAndGet();
                    condition.signalAll();
                } catch (InterruptedException ignored) {}
                finally {
                    lock.unlock();
                }
            }, "t0");
        }

        for (int i = 0; i < threadCount; ++i) {
            threads[i].start();
        }

        try {
            for (int i = 0; i < threadCount; ++i) {
                threads[i].join();
            }
        } catch (InterruptedException ignored) {}


        final String[] cubeShow = new String[1];
        assertDoesNotThrow(() -> { cubeShow[0] = cube.show(); });
        assertTrue(countCubeColor(cubeShow[0], size));

        assertEquals(afterCounter.get(), beforeCounter.get());

        var rapport = insideRapport.stream().collect(Collectors.groupingBy( Function.identity(), Collectors.counting()));
        assertTrue(rapport.containsKey(0) && rapport.size() == 1);

        assertFalse(error[0]);

        assertEquals(0, logger.size() % 2);

        Move[] moveLogs = logger.toArray(new Move[0]);
        for (int i = 0; i < logger.size(); i += 2) {
            assertEquals(moveLogs[i], moveLogs[i + 1]);
        }
        for (int i = logger.size() - 1; i > 0 ; --i) {
            assertTrue(moveLogs[i].beforeOfEqualByLayer(moveLogs[i - 1]));
        }
    }

    @Test
    @DisplayName("Small Interrupt Test")
    public void testSmallInterrupt() {
        final int size = 1000;
        Cube cube = generateSimpleCube(size);
        final AtomicBoolean good = new AtomicBoolean(true);

        Thread t1 = new Thread(() -> {
            try {
                int t = size;
                while (t --> 0) {
                    cube.rotate(0, 0);
                }
                good.set(false);
            } catch (InterruptedException ignored) {}
        } ,"t1");

        Thread t2 = new Thread(() -> {
            try {
                int t = size;
                while (t --> 0) {
                    cube.rotate(0, 0);
                }
            } catch (InterruptedException e) {
                good.set(false);
            }
        } ,"t2");

        t1.start();
        t2.start();

        try {
            Thread.sleep(1000);
            t1.interrupt(); // The cube is huge, so the thread should not be able to finish the task.

            t1.join();
            t2.join();
        } catch (InterruptedException ignored) {}

        assertTrue(good.get());
    }
}
