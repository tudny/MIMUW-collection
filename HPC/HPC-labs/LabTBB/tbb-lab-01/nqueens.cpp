// Solves the n-quees puzzle on an n x x checkerboard.
//
// This sequential implementation is to be extended with TBB to get a
// parallel implementation.
//
// HPC course, MIM UW
// Krzysztof Rzadca, LGPL

#include "tbb/tbb.h"
#include <iostream>
#include <list>
#include <vector>
#include <cmath>


#define THREAD_DEPTH 4


// Indexed by column. Value is the row where the queen is placed,
// or -1 if no queen.
typedef std::vector<int> Board;


void pretty_print(const Board& board) {
    for (int row = 0; row < (int) board.size(); row++) {
        for (const auto& loc : board) {
            if (loc == row)
                std::cout << "*";
            else
                std::cout << " ";
        }
        std::cout << std::endl;
    }
}


// Checks the location of queen in column 'col' against queens in cols [0, col).
bool check_col(Board& board, int col_prop) {
    int row_prop = board[col_prop];
    int col_queen = 0;
    for (auto i = board.begin();
         (i != board.end()) && (col_queen < col_prop);
         ++i, ++col_queen) {
        int row_queen = *i;
        if (row_queen == row_prop) {
            return false;
        }
        if (abs(row_prop - row_queen) == col_prop - col_queen) {
            return false;
        }
    }
    return true;
}


void initialize(Board& board, int size) {
    board.reserve(size);
    for (int col = 0; col < size; ++col)
        board.push_back(-1);
}


void initialize_parallel(Board& board, int size) {
    board.clear();
    board.reserve(size);
}


// Solves starting from a partially-filled board (up to column col).
void recursive_solve(Board& partial_board, int col, std::list<Board>& solutions) {
    // std::cout << "rec solve col " << col << std::endl;
    // pretty_print(b_partial);

    int b_size = partial_board.size();
    if (col == b_size) {
        solutions.push_back(partial_board);
    }
    else {
        for (int tested_row = 0; tested_row < b_size; tested_row++) {
            partial_board[col] = tested_row;
            if (check_col(partial_board, col))
                recursive_solve(partial_board, col+1, solutions);
        }
    }
}


void parallel_solve(Board& partial_board, int b_size, std::vector<Board> &result, std::atomic<int>& solution_count) {
    tbb::enumerable_thread_specific<std::vector<Board>> solutions;
    tbb::parallel_for_each(&partial_board, &partial_board + 1, [&](const Board &board, tbb::feeder<Board>& feeder) {
        int col = board.size();
        if (col == b_size) {
            solutions.local().push_back(board);
            ++solution_count;
            return;
        }

        for (int tested_row = 0; tested_row < b_size; tested_row++) {
            Board new_board = board;
            new_board.push_back(tested_row);
            if (check_col(new_board, col)) {
                feeder.add(new_board);
            }
        }
    });
    result.reserve(solution_count);
    for (auto &sub_result : solutions) {
        result.insert(result.end(), sub_result.begin(), sub_result.end());
    }
}


int main() {
    const int board_size = 13;
    Board board{};
    initialize(board, board_size);
    std::list<Board> solutions{};

    tbb::tick_count seq_start_time = tbb::tick_count::now();
    recursive_solve(board, 0, solutions);
    tbb::tick_count seq_end_time = tbb::tick_count::now();
    double seq_time = (seq_end_time-seq_start_time).seconds();
    std::cout << "seq time: " << seq_time << "[s]" <<std::endl;

    std::cout << "solution count: "<<solutions.size()<<std::endl;
    // for (const auto& sol : solutions) {
    //     pretty_print(sol);
    //     std::cout << std::endl;
    // }

    board.clear();
    initialize_parallel(board, board_size);
    std::vector<Board> solutions_queue_list{};
    std::atomic<int> solution_count{0};
    tbb::tick_count par_start_time = tbb::tick_count::now();
    parallel_solve(board, board_size, solutions_queue_list, solution_count);
    tbb::tick_count par_end_time = tbb::tick_count::now();
    double par_time = (par_end_time-par_start_time).seconds();
    std::cout << "par time: " << par_time << "[s]" <<std::endl;

    std::cout << "solution count: "<< solutions.size() << std::endl;
    std::cout << "solution count: "<< solution_count << std::endl;
    // for (const auto& sol : solutions) {
    //     pretty_print(sol);
    //     std::cout << std::endl;
    // }

    // Compare the same
    std::sort(solutions_queue_list.begin(), solutions_queue_list.end());

    if (solutions.size() != solutions_queue_list.size()) {
        std::cout << "Different solution count" << std::endl;
    } else {
        auto it1 = solutions.begin();
        auto it2 = solutions_queue_list.begin();
        while (it1 != solutions.end()) {
            if (*it1 != *it2) {
                std::cout << "Different solution" << std::endl;
                break;
            }
            it1++;
            it2++;
        }
    }

    return 0;
}
