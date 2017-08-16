#!/usr/bin/python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import glob
import os
import sys


def main():
    parser = OptionParser()
    usage = "Usage: %prog [options] arg1 arg2"
    parser.add_option("-f", "--folder", dest="folder",
                      help="folder of read files")
    parser.add_option("-d", "--design", dest="design",
                      help="output design file")
    parser.add_option("-p", "--paired", default=False, dest="paired",
                      action="store_true",
                      help="tag to indicate the data are paired end")
    (options, args) = parser.parse_args()
    option_dict = vars(options)
    folder = option_dict['folder']
    design = option_dict['design']
    paired = option_dict['paired']
    if folder is None or design is None:
        parser.error("options -f and -d should be set")
    files1 = glob.glob(os.path.join(folder, '*.fastq.gz'))
    files2 = glob.glob(os.path.join(folder, '*.fq.gz'))
    files1.extend(files2)
    files = sorted(files1)
    with open(design, 'w') as file_design:
        output_list =[]
        if paired  :
            sample_num = int(len(files)/2)
            if sample_num != len(files)/2:
                sys.exit("Odd numbers of read files in paired design\n")
            else:
                head = ("sample", "left", "right")
                file_design.write("\t".join(head) + "\n")
                for count in range(sample_num):
                    file1 = files[count*2]
                    file2 = files[count*2+1]
                    dir, file_short = os.path.split(file1)
                    output_list.append([file_short,file1 , file2])
        else:
            head = ("sample", "left")
            file_design.write("\t".join(head) + "\n")
            sample_num = len(files)
            for count in range(sample_num):
                file1 = files[count]
                dir, file_short  = os.path.split(file1)
                output_list.append([file_short, file1])
        # change sample name
        for line_list in output_list:
            fq1 = line_list[0]
            checking = ("_R1_001.fastq.gz","_R1_001.fq.gz","_1.fastq.gz", "_1.fq.gz",
                        "_R1.fastq.gz","_R1.fq.gz",".fastq.gz",".fq.gz")
            for suffix in checking:
                if fq1.endswith(suffix):
                    sample = fq1[0:((-1)*len(suffix))]
                    line_list[0] = sample
                    break
            file_design.write("\t".join(line_list)+"\n")
#######################################################
#  Main function
#######################################################

if __name__ == '__main__':
    main()



