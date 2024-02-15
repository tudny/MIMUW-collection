#!/usr/bin/env python3

import os
import sys

WHAT_TO_CHANGE = """
%[CHANGETHIS]
program :-
    % Ensure list library
    ensure_loaded(library(lists)),
%[ENDOFCHANGE]
"""

CHANGE_TO = """
:- ensure_loaded(library(lists)).

user:runtime_entry(start) :-
"""


def main():
    # Open file given as argv
    filename = sys.argv[1]
    with open(filename, "r") as f:
        lines = f.readlines()
        lines = ''.join(lines)

    # Replace WHAT_TO_CHANGE with CHANGE_TO
    lines = lines.replace(WHAT_TO_CHANGE, CHANGE_TO)
    with open(f'students/{filename}', "w") as f:
        f.write(lines)


if __name__ == '__main__':
    main()
    print('Done!')
