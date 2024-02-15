from abc import ABC

import numpy as np
import matplotlib.pyplot as plt


def distance(s1: str, s2: str) -> int:
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


class ChainType:
    def __init__(self):
        pass

    def distance(self, s1: str, s2: str) -> int:
        raise NotImplementedError

    def get_matrix(self) -> dict:
        raise NotImplementedError


class JC69(ChainType, ABC):
    def __init__(self, miu: float):
        super().__init__()
        self.miu = miu

    def distance(self, s1: str, s2: str) -> int:
        p = distance(s1, s2) / len(s1)
        return -3 / 4 * np.log(1 - 4 / 3 * p)

    def get_matrix(self) -> dict:
        miu = self.miu
        acgt = "ACGT"
        xy = 1 / len(acgt) * miu
        xx = 1 - (len(acgt) - 1) * xy
        return {c: {c2: (xx if c == c2 else xy) for c2 in acgt} for c in acgt}


class MarkovChain:
    def __init__(self, start: str, model: ChainType):
        self.state = start
        self.chain = model.get_matrix()

    def _map_single(self, c: str) -> str:
        return np.random.choice(
            list(self.chain[c].keys()), p=list(self.chain[c].values())
        )

    def step(self):
        self.state = "".join(self._map_single(c) for c in self.state)

    def simulate(self, n: int, callback: callable = None):
        for _ in range(n):
            self.step()
            if callback is not None:
                callback(self.state)


def plot(model_distances, simple_distances, title, file_name):
    fig, ax1 = plt.subplots()
    color1 = "tab:red"
    ax1.set_xlabel("Model steps")
    ax1.set_ylabel("Model distance", color=color1)
    ax1.plot(model_distances, color=color1)
    ax1.tick_params(axis="y", labelcolor=color1)

    ax2 = ax1.twinx()

    color2 = "tab:blue"
    ax2.set_ylabel("Simple distance", color=color2)
    ax2.plot(simple_distances, color=color2)
    ax2.tick_params(axis="y", labelcolor=color2)

    fig.tight_layout()
    plt.title(title)
    plt.savefig(f"data/{file_name}.svg")
    plt.show()


def main():
    start_len = 1000
    steps = 1000
    mius = [10e-10, 10e-5, 10e-3, 10e-2]
    # mius = [1e-3]

    start = "".join(np.random.choice(list("ACGT"), start_len))

    import os

    if not os.path.exists("data"):
        os.makedirs("data")

    for miu in mius:
        jc_model = JC69(miu)
        mc = MarkovChain(start, jc_model)
        distances = []
        simple_distances = []

        def append_distance(s):
            distances.append(jc_model.distance(start, s))
            simple_distances.append(distance(start, s))

        mc.simulate(steps, append_distance)
        plot(distances, simple_distances, f"JC69 model with {miu=}", f"jc69-{miu}")
        distances.clear()


if __name__ == "__main__":
    main()
