import os
import sys

from Parser import Parser


def check_ars_validity(args):
    """
    Check that the args are in the format of :
    <input_file_1.out> <input_file_2.out> ... <input_file_n.out> <output_file.xyz>
    :param args:
    :return: True if the args are in the correct format and false otherwise
    """
    if len(args) == 0 or len(args) == 1:
        print('Please call the program with an input file/s name/s and an output file name with xyz extension')
        return False
    if not args[-1].endswith('.xyz'):
        print('The output file should end with .xyz extension.')
        return False
    for file in args[:-1]:
        pass
    return True


def parse_files(input_files_names, parser):
    """
    Recursively parse each file in the given directory or files
    :param input_files_names: name of directories and files to parse
    :param parser: the current parser to parse with
    """
    for input_file_name in input_files_names:
        path = os.path.abspath(input_file_name)
        if os.path.isdir(path):
            # recursively parse each file in the directory
            parse_files([os.path.join(path, file_name) for file_name in os.listdir(input_file_name)], parser)
        else:
            input_file = open(input_file_name, 'r')
            parser.set_new_file(input_file)
            parser.parse_file()
            input_file.close()


def main(args):
    """
    Go over all the given files and parse them into one xyz file.
    :param args: a list on file names to parse
    """
    if not check_ars_validity(args):
        return
    input_files_names, output_file_name = args[:-1], args[-1]
    output_file = open(output_file_name, 'w')
    # the origin_list_file will have a list of names of the files that were parsed
    origin_list_file = open(output_file_name.replace('xyz', 'txt'), 'w')
    parser = Parser(output_file, origin_list_file)
    parse_files(input_files_names, parser)
    output_file.close()
    origin_list_file.close()


if __name__ == '__main__':
    main(sys.argv[1:])
