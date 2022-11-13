#include <utility>
#include <deque>
#include <future>
#include <thread>

#include "teams.hpp"
#include "contest.hpp"
#include "collatz.hpp"

class ConcurrentQueue {
    std::queue<size_t> queue;
    std::mutex mutex;

public:
    void push(const size_t &element) {
        std::lock_guard<std::mutex> lock(mutex);
        queue.push(element);
    }

    size_t front() {
        std::lock_guard<std::mutex> lock(mutex);
        return queue.front();
    }

    __attribute__((unused)) std::queue<std::size_t>::size_type size() {
        std::lock_guard<std::mutex> lock(mutex);
        return queue.size();
    }

    bool empty() {
        std::lock_guard<std::mutex> lock(mutex);
        return queue.empty();
    }

    void pop() {
        std::lock_guard<std::mutex> lock(mutex);
        queue.pop();
    }
};

ContestResult TeamNewThreads::runContestImpl(ContestInput const &contestInput) {
    ContestResult r;

    const auto sizeOfInput = contestInput.size();
    const auto workersLimit = getSize();
    std::atomic_int workersRunningCounter{0};
    std::condition_variable tooManyThreads;
    std::mutex mutex;

    std::vector<std::promise<uint64_t>> resultPromises(sizeOfInput);
    std::vector<std::future<uint64_t>> resultFutures(sizeOfInput);
    std::vector<std::thread> workers(sizeOfInput);
    ConcurrentQueue threadToAwait;
    std::vector<bool> alreadyJoined(sizeOfInput, false);

    for (size_t id = 0; id < sizeOfInput; ++id) {
        resultFutures[id] = resultPromises[id].get_future();
    }

    auto workerJob = [&workersRunningCounter, &tooManyThreads, &threadToAwait]
            (size_t id, InfInt const &input, std::promise<uint64_t> &resultPromise) {

        auto result = calcCollatz(input);
        resultPromise.set_value(result);

        --workersRunningCounter;

        threadToAwait.push(id);

        tooManyThreads.notify_all();
    };

    std::unique_lock<std::mutex> mutexLock(mutex);
    for (size_t id = 0; id < sizeOfInput; ++id) {
        tooManyThreads.wait(mutexLock, [&workersRunningCounter, &workersLimit](){ return workersRunningCounter < workersLimit; });

        if (!threadToAwait.empty()) {
            workers[threadToAwait.front()].join();
            alreadyJoined[threadToAwait.front()] = true;
            threadToAwait.pop();
        }

        ++workersRunningCounter;
        workers[id] = createThread(workerJob, id, std::ref(contestInput[id]), std::ref(resultPromises[id]));
    }

    r = cxxpool::get(resultFutures.begin(), resultFutures.end());

    for (size_t id = 0; id < sizeOfInput; ++id) {
        if (!alreadyJoined[id])
            workers[id].join();
    }

    return r;
}

ContestResult TeamConstThreads::runContestImpl(ContestInput const &contestInput) {
    ContestResult r;

    const auto numberOfThreads = getSize();
    const auto sizeOfInput = contestInput.size();

    std::vector<std::thread> workers(numberOfThreads);
    std::vector<std::promise<uint64_t>> resultPromises(sizeOfInput);
    std::vector<std::future<uint64_t>> result_futures(sizeOfInput);

    for (size_t id = 0; id < sizeOfInput; ++id) {
        result_futures[id] = resultPromises[id].get_future();
    }

    auto worker_job = [&contestInput, &resultPromises, &sizeOfInput, &numberOfThreads](size_t id) {
        for (size_t index = id; index < sizeOfInput; index += numberOfThreads) {
            resultPromises[index].set_value(calcCollatz(contestInput[index]));
        }
    };

    for (size_t id = 0; id < numberOfThreads; ++id) {
        workers[id] = createThread(worker_job, id);
    }

    r = cxxpool::get(result_futures.begin(), result_futures.end());

    for (auto &worker : workers) {
        worker.join();
    }

    return r;
}

ContestResult TeamPool::runContest(ContestInput const &contestInput) {
    ContestResult r;
    //TODO
    return r;
}

ContestResult TeamNewProcesses::runContest(ContestInput const &contestInput) {
    ContestResult r;
    //TODO
    return r;
}

ContestResult TeamConstProcesses::runContest(ContestInput const &contestInput) {
    ContestResult r;
    //TODO
    return r;
}

ContestResult TeamAsync::runContest(ContestInput const &contestInput) {
    ContestResult r;
    //TODO
    return r;
}
