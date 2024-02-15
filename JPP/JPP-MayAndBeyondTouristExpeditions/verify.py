#!/usr/bin/env python3

import os
import re


def error(message):
    print(message)
    exit(1)


def get_filename_from_args():
    import sys
    if len(sys.argv) < 2:
        return None
    return sys.argv[1]


def get_filename_from_input():
    filename = input("Enter filename: ")
    return filename


filename_regex = r'([a-zA-Z])([a-zA-Z])[0-9]{6}\.pl'

def check_file_name_is_index_with_pl_extension(filename):
    match = re.fullmatch(filename_regex, filename)
    if match:
        return match.group(1), match.group(2)
    return None


def check_if_file_exists(filename):
    if not os.path.isfile(filename):
        print("File does not exist")
        return False
    return True


def check_if_all_lines_in_file_are_less_than_80_characters(filename):
    with open(filename, "r") as f:
        for id, line in enumerate(f):
            if len(line) - 1 > 80:     # We don't count new line char.
                print("Line {} is too long".format(id + 1))
                return False
    return True


surename_regex = r'\s*%\s*([A-Za-z]*-?[A-Za-z]+)\s+([A-Za-z]*-?[A-Za-z]+)\s*\n'

def check_if_first_line_in_file_contains_a_comment_with_surename(filename, name_letter, surename_letter):
    with open(filename, "r") as f:
        first_line = f.readline()
        if not first_line.startswith("%"):
            print("First line does not contain a comment")
            return False
        if not (match := re.fullmatch(surename_regex, first_line)):
            print("First line does not contain a comment with a surename")
            return False
        if not match.group(1).lower().startswith(name_letter.lower()) or\
              not match.group(2).lower().startswith(surename_letter.lower()):
            print("First line does not contain a comment with a surename starting with {}".format(surename_letter))
            return False
    return True


def main():
    filename = get_filename_from_args() or get_filename_from_input()
    if not check_if_file_exists(filename):
        error("File does not exist")
    if not (names := check_file_name_is_index_with_pl_extension(filename)):
        error("File name is not index with pl extension")
    (name_letter, surename_letter) = names
    if not check_if_all_lines_in_file_are_less_than_80_characters(filename):
        error("Line is too long")
    if not check_if_first_line_in_file_contains_a_comment_with_surename(filename, name_letter, surename_letter):
        error("First line does not contain a comment with a surename")
    print("File is OK")
    

if __name__ == "__main__":
    main()
