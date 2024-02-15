def sqGen(n):
    sq = 0
    for k in range(1, n + 1):
        odd = 2 * k - 1
        sq += odd
        yield sq


for i in sqGen(10):
    print(i)
