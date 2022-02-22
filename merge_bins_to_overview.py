#!/usr/bin/env python
import sys

# define input and output files
input_txt_file = sys.argv[1]
output_txt_file = input_txt_file.replace(".txt", "_bin_merged.txt")
input_bin_file = sys.argv[2]

# create a dictionary of scaffold to bin relation
bin_dict = dict()
with open(input_bin_file, 'r') as binfile:
    for line in binfile:
        line_tuple = tuple(filter(None, line.rstrip().split("\t")))
        if line_tuple[0] in bin_dict:
            print("error multiple bins for " + line_tuple[0])
        else:
            bin_dict[line_tuple[0]] = line_tuple[1]

# cross reference scaffold name to above dictionary and create a new file
with open(output_txt_file, "w+") as outprint:
    with open(input_txt_file, 'r') as txtfile:
        header = txtfile.readline().rstrip() + "\t" + "bin\n"
        outprint.write(header)
        for line in txtfile:
            line_tuple = tuple(filter(None, line.rstrip().split("\t")))
            bine = ""
            ine = line.rstrip() + "\t"
            if line_tuple[0] in bin_dict:
                bine = bin_dict[line_tuple[0]]
                line = line.rstrip() + "\t" + bine + "\n"
            else:
                line = line.rstrip() + "\t\n"
            outprint.write(line)
