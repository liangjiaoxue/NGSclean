#!/usr/bin/python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import glob
import os
import gzip
import sys


class ReadRecord:
    def __init__(self,thread):
        self.paired = False
        self.thread = thread
        self.samples = ()
        self.sample2len = dict()
        self.sample2phred = dict()

    def read_design(self, design_file):
        with open(design_file, 'r') as design_f:
            lines = [line.rstrip().split("\t") for line in design_f]
            if len(lines[0]) == 3:
                self.paired = True
            elif len(lines[0]) == 2:
                self.paired = False
            else:
                print("Incorrect design format, please double check the file")
            if self.paired :
                if lines[0][2] != 'right':
                    sys.exit("Please check the head of design files")
            self.samples = lines[1:]

    def check_fastq(self):
        for sample_record in self.samples:
            first_fqgz = sample_record[1]
            sample_name = sample_record[0]
            sample_name.replace(' ', '_')
            # get the phred format
            len_max = 0
            qual_min = 1000
            qual_max = 0
            with gzip.open(first_fqgz, "rt") as fqgzIN:
                count = 0
                limit = 200
                qual_str_all = ''
                for line in fqgzIN:
                    count += 1
                    if count > limit:
                        break
                    if count % 4 == 0:
                        qual_str = line.rstrip()
                        qual_str_all += qual_str
                        if len_max < len(qual_str):
                            len_max = len(qual_str)
                for qual_char in qual_str_all:
                    qual = ord(qual_char)
                    if qual_min > qual:
                        qual_min = qual
                    if qual_max < qual:
                        qual_max = qual
            if qual_min < 33 or qual_max > 105:
                sys.exit("Quality values corrupt. found [$min; $max] where [33; 104] was expected \n")
            if qual_min >= 59 and qual_max <= 110:
                qual_out = '64'
            elif qual_min >= 33 and qual_max <= 74:
                qual_out = '33'
            else:
                sys.exit("May be new fastq format \n")
            # output to class
            self.sample2len[sample_name] = len_max
            self.sample2phred[sample_name] = qual_out

def main():
    ########################################################################
    # Options
    parser = OptionParser()
    usage = "Usage: %prog [options] arg1 arg2"
    parser.add_option("-d", "--design", dest="design",
                      help="Input design file")
    parser.add_option("-s", "--singleton", dest="singleton", default="keep",
                      help="How to treat singleton reads: keep(default), merge with paired, or discard")
    parser.add_option("-t", "--threadNum", dest="thread",
                      help="Numbers of thread")
    parser.add_option("-o", "--outputDir", dest="output",
                      help="Change the workingDir for read output")
    parser.add_option("-q", "--queue", dest="queue",
                      help="Queue for batch jobs, inter (run interactively) or queue (default)")
    parser.add_option("-n", "--node", dest="node",
                      help="Node of queue")
    parser.add_option("-z", "--run_trimmomatic", dest="trimmomatic",
                      help="path to trimmomatic.jar")
    parser.add_option("-x", "--load_trimmo_module", dest="loadtrimmo",
                      help="module for trimmomatic")
    parser.add_option("-a", "--adaptor", dest="adaptor",
                      help="path to trimmomatic adaptor")
    parser.add_option("-y", "--run_star", dest="star",
                      help="path to STAR")
    parser.add_option("-w", "--load_star_module", dest="loadstar",
                      help="path to STAR")
    #parser.add_option("-r", "--rrna", dest="rrna",
                 #     help="Reference of rRNA")

    (options, args) = parser.parse_args()
    # parser options
    option_dict = vars(options)
    design = option_dict['design']
    singleton  = option_dict['singleton']
    thread = option_dict['thread']
    #rrna = option_dict['rrna']
    output = option_dict['output']
    queue  = option_dict['queue']
    node   = option_dict['node']
    adaptor= option_dict['adaptor']
    trimmomatic = option_dict['trimmomatic']
    loadtrimmo = option_dict['loadtrimmo']
    loadstar  = option_dict['loadstar']
    star = option_dict['star']

    #
    if thread is None:
        thread = 1
    if queue is None:
        queue = 'batch'
    if node is None:
        node = 'HIGHMEM'
    if trimmomatic is None or star is None:
        sys.exit("Please provide the path to trimmomatic.jar and STAR")
    if adaptor is None:
        sys.exit("Please provide the path to trimmomatic adaptor")
    #if rrna is None :
        #sys.exit("Please provide the path to STAR reference of rRNA sequences")
    if output is None:
        output = os.getcwd()
    python_dir_path = os.path.dirname(os.path.realpath(__file__))
    #print(python_dir_path +"path to python code")
    rrna = python_dir_path + '/rRNA_ref'
    ###########################################################################
    # read design file
    readrecord = ReadRecord(thread)
    readrecord.read_design(design)
    readrecord.check_fastq()

    # prepare shell script
    dir, file_design_short = os.path.split(design)
    master_file = 'Run_'+file_design_short[0:-4]+'.sh'
    working_dir = output

    master_out = open(master_file, 'w')
    master_out.write("#!/bin/sh\n")
    master_out.write("cd " + working_dir + "\n")

    #############################################################################
    ## Each sample
    count = 0
    for sample_record in readrecord.samples:
        ## get children shell name
        sample_name = sample_record[0]
        sample_name.replace(' ', '_')
        first_fastq = sample_record[1]
        first_fastq.replace(' ', '')
        dir, first_fastq_short = os.path.split(first_fastq)
        prefix = sample_name
        #prefix = first_fastq_short
        #if 'fastq.gz' in first_fastq_short:
        #    prefix = first_fastq_short[0:-9]
        #elif 'fq.gz' in first_fastq_short:
        #    prefix = first_fastq_short[0:-6]

        count += 1
        file_shell_str = 'r' + str(count) + '_' + prefix + '.sh'
        if queue == 'inter':
            master_out.write("chmod 750 " + file_shell_str + "\n")
            master_out.write("./" + file_shell_str + "\n")
        else:
            master_out.write("qsub " + file_shell_str + "\n")

        ##############################
        ## write childredn shell
        write_shell_head(file_shell_str, queue, thread, node, working_dir)

        ##############################
        ## write trmmomatic shell
        all_fastq_file = " ".join(sample_record[1:])
        write_shell_trimmo(file_shell_str, loadtrimmo, trimmomatic, prefix, all_fastq_file,readrecord,sample_name,adaptor)

        ###############################
        ### STAR module
        file_shell = open(file_shell_str, 'a')
        if loadstar is None:
            print("STAR module is not loaded")
        else:
            file_shell.write("module load " + loadstar + "\n")
        file_shell.close()
        # STAR
        mis_cut = int(readrecord.sample2len[sample_name] * 0.1)
        # mismatch cutoff for rRNA
        # checking paried
        star_out_prefix = prefix + '_STAR_'
        if readrecord.paired:
            # Paired end
            r1_pair = prefix + "_trimP_1.fq.gz"
            r2_pair = prefix + "_trimP_2.fq.gz"
            star_input = r1_pair + ' ' + r2_pair
            write_shell_star(file_shell_str, star, star_input, star_out_prefix, thread, rrna, mis_cut)
            if singleton != 'discard':
                r1_unpair = prefix + "_trimS_1.fq.gz"
                r2_unpair = prefix + "_trimS_2.fq.gz"
                # single 1
                star_input = r1_unpair
                star_out_prefix = prefix + '_STAR_S1_'
                write_shell_star(file_shell_str, star, star_input, star_out_prefix, thread, rrna, mis_cut)
                # single 2
                star_input = r2_unpair
                star_out_prefix = prefix + '_STAR_S2_'
                write_shell_star(file_shell_str, star, star_input, star_out_prefix, thread, rrna, mis_cut)
        else:
            # single end
            pe_tag = 'SE'
            star_input = prefix + '_trim_1.fq.gz'
            write_shell_star(file_shell_str,star, star_input, star_out_prefix, thread, rrna, mis_cut)


        ###################################################
        ### Summary data
        final_log = prefix + "_Final.log.txt"
        trim_log = prefix + '_trim.log'
        star_main_prefix = prefix + '_STAR_'
        write_shell_summary(file_shell_str,final_log,trim_log,star_main_prefix)

        ###################################################
        ###  prepare fastq output
        paired_boolean = readrecord.paired
        write_shell_fastq(file_shell_str, paired_boolean, prefix,singleton)

        ##################################################
        ### zip and clean
        write_shell_clean(file_shell_str,star_main_prefix,prefix, paired_boolean,singleton)
# End of main function
################################################################################################################


def  write_shell_clean(file_shell_str,star_main_prefix,prefix, paired_boolean,singleton):
    file_shell = open(file_shell_str, 'a')
    remove_list = []
    output_list = []
    star_temp_files = ("SJ.out.tab", "Log.progress.out", "Log.out", "Log.final.out", "Aligned.out.sam")
    if not paired_boolean:
        # single end
        # STAR
        for suffix in star_temp_files:
            remove_list.append(star_main_prefix + suffix)
        # Trim
        remove_list.append(prefix + '_trim_1.fq.gz')
        remove_list.append(prefix + '_trim.log')
        output_list.append(prefix + '_clean_1.fq')
    else:
        # paired
        ## Main STAR
        for suffix in star_temp_files:
            remove_list.append(star_main_prefix + suffix)
        ## single end STAR
        if singleton != 'discard':
            star_out_prefix = prefix + '_STAR_S1_'
            for suffix in star_temp_files:
                remove_list.append(star_out_prefix + suffix)
            star_out_prefix = prefix + '_STAR_S2_'
            for suffix in star_temp_files:
                remove_list.append(star_out_prefix + suffix)

        if singleton == 'merge':
            fastq_temp_list = (
                "_clean_S1.fq", "_clean_S2.fq", "_clean_E1.fq", "_clean_E2.fq", "_STAR_S2_Unmapped.out.mate1")
            for suffix in fastq_temp_list:
                remove_list.append(prefix + suffix)
        elif singleton == 'keep':
            remove_list.append(prefix + '_STAR_S2_Unmapped.out.mate1')
            output_list.append(prefix + '_clean_S1.fq')
            output_list.append(prefix + '_clean_S2.fq')
        elif singleton == 'discard':
            pass
        else:
            sys.exit("Wrong method to handle singleton reads\n")
        # Trim
        trim_temp_files = ("_STAR_Unmapped.out.mate1","_STAR_Unmapped.out.mate2",
            "_trimP_1.fq.gz", "_trimS_1.fq.gz", "_trimP_2.fq.gz", "_trimS_2.fq.gz", "_trim.log")
        for suffix in trim_temp_files:
            remove_list.append(prefix + suffix)
        # Output
        output_list.append(prefix + '_clean_1.fq')
        output_list.append(prefix + '_clean_2.fq')
    ## Write shell
    for remove_f in remove_list:
        file_shell.write("rm " + remove_f + "\n")
    for gfile in output_list:
        file_shell.write("gzip -f " + gfile + "\n")
    file_shell.close()


def write_shell_head(file_shell_str, queue,thread,  node, working_dir):
    # children shell  head
    file_shell = open(file_shell_str, 'w')
    file_shell.write("#!/bin/sh\n" +
                "#PBS -N j-clean" + file_shell_str[0:-3] + "\n" +
                "#PBS -q " + queue + "\n" +
                "#PBS -l nodes=1:ppn=" + str(thread) + ':' + node + "\n" +
                "#PBS -l walltime=48:00:00\n" +
                "#PBS -l mem=20gb\n" +
                "cd " + working_dir + "\n\n")
    file_shell.close()

def write_shell_trimmo(file_shell_str, loadtrimmo, trimmomatic, prefix, all_fastq_file,readrecord,sample_name,adaptor):
    file_shell = open(file_shell_str, 'a')
    paired_boolean = readrecord.paired
    phred_tag = ''
    if readrecord.sample2phred[sample_name] == '64':
        phred_tag = '-phred64'
    else:
        phred_tag = '-phred33'
    paired_boolean = readrecord.paired
    thread = readrecord.thread
    ### write
    if loadtrimmo is None:
        print("trimmomatic module is not loaded")
    else:
        file_shell.write("module load " + loadtrimmo + " \n")
    pe_tag = 'SE'
    trim_out = ''
    trim_log = prefix + '_trim.log'
    if paired_boolean:
        # Paired end
        pe_tag = 'PE'
        r1_pair = prefix + "_trimP_1.fq.gz"
        r1_unpair = prefix + "_trimS_1.fq.gz"
        r2_pair = prefix + "_trimP_2.fq.gz"
        r2_unpair = prefix + "_trimS_2.fq.gz"
        trim_out_list = (r1_pair, r1_unpair, r2_pair, r2_unpair)
        trim_out = " ".join(trim_out_list)
    else:        # single end
        pe_tag = 'SE'
        trim_out = prefix + '_trim_1.fq.gz'
    file_shell.write("time java -jar " + trimmomatic + ' ' + pe_tag + " -threads " + str(thread) + " " +
                     phred_tag + " \\\n" +
                     " " + all_fastq_file + " \\\n" +
                     " " + trim_out + " \\\n" +
                     " ILLUMINACLIP:" + adaptor + ':2:30:10' +
                     ' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 &>' + trim_log + "\n\n")
    file_shell.close()


def write_shell_star(file_shell_str, star, star_input, star_out_prefix, thread, rrna, mis_cut):
    with open(file_shell_str, "a") as shell_out:
        shell_out.write(star + ' --runThreadN ' + str(thread) + ' --genomeDir ' + rrna + " \\\n" +
                        " --readFilesIn " + star_input + " \\\n" +
                        " --readFilesCommand gunzip -c --outReadsUnmapped Fastx " + " \\\n" +
                        " --outFileNamePrefix " + star_out_prefix + "  --outFilterMismatchNmax " + str(
            mis_cut) + " \\\n" +
                        " --alignMatesGapMax 20 --alignIntronMax 20 \n")
        shell_out.write("\n")


def write_shell_summary(file_shell_str, final_log, trim_log, star_prefix):
    file_shell = open(file_shell_str, 'a')
    file_shell.write('printf "Trimming : " >' + final_log + "\n")
    file_shell.write('grep "^Input Read"  ' + trim_log + '>>' + final_log + "\n")
    file_shell.write('printf "Mapping rRNA : " >>' + final_log + "\n")
    file_shell.write(
        "grep \"Number of input read\"  " + star_prefix + "Log.final.out >>" + final_log + "\n")
    file_shell.write('printf "Mapping rRNA : " >>' + final_log + "\n")
    file_shell.write(
        "grep \"Uniquely mapped reads number\"  " + star_prefix + "Log.final.out >>" + final_log + "\n")
    file_shell.write('printf "Mapping rRNA : " >>' + final_log + "\n")
    file_shell.write(
        "grep \"Number of reads mapped to multiple loci\"  " + star_prefix + "Log.final.out >>" + final_log + "\n")
    file_shell.write('printf "HPT : " >>' + final_log + "\n")
    file_shell.write(
        'grep "MARKER" ' + star_prefix + 'Aligned.out.sam |grep -v "^@"|grep "MARKER_HPT"|cut -f 1|sort|uniq|wc -l' +
        '>>' + final_log + "\n")
    file_shell.write('printf "NPT : " >>' + final_log + "\n")
    file_shell.write(
        'grep "MARKER" ' + star_prefix + 'Aligned.out.sam |grep -v "^@"|grep "MARKER_NPT"|cut -f 1|sort|uniq|wc -l' +
        '>>' + final_log + "\n\n")
    file_shell.close()


def write_shell_fastq(file_shell_str, paired_boolean, prefix, singleton):
    star_prefix = prefix + '_STAR_'
    file_shell = open(file_shell_str, 'a')
    if not paired_boolean:
        # single send
        final_out1 = prefix + "_clean_1.fq"
        file_shell.write("mv " + star_prefix + "Unmapped.out.mate1 " + final_out1 + "\n")
    else:
        # Paired end
        final_out1 = prefix + "_clean_1.fq"
        final_out2 = prefix + "_clean_2.fq"
        #file_shell.write("mv " + star_prefix + "Unmapped.out.mate1 " + final_out1 + "\n")
        #file_shell.write("mv " + star_prefix + "Unmapped.out.mate2 " + final_out2 + "\n")
        file_shell.write(''' sed '1~4 s|\s00$||g' < ''' + star_prefix + "Unmapped.out.mate1" + ' >' + final_out1 + "\n")
        file_shell.write(''' sed '1~4 s|\s00$||g' < ''' + star_prefix + "Unmapped.out.mate2" + ' >' + final_out2 + "\n")
        #####################################
        read1_ori = prefix + '_STAR_S1_' + 'Unmapped.out.mate1'
        read2_ori = prefix + '_STAR_S2_' + 'Unmapped.out.mate1'
        if singleton != 'discard':
            read1_rename = prefix + "_clean_S1.fq"
            file_shell.write("mv " + read1_ori + " " + read1_rename + "\n")
            read2_rename = prefix + "_clean_S2.fq"
            file_shell.write(''' sed '1~4 s|/1$|/2|g' < ''' + read2_ori + ' >' + read2_rename + "\n")
            if singleton == 'merge':
                # regenerate records
                fake_fq2_from_fq1 = prefix + "_clean_E2.fq"
                file_shell.write(''' sed '1~4 s|/1$|/2|g' < ''' + read1_rename + ' | ' + " \\\n" +
                                 ''' awk 'NR%4==1{print}NR%4==2{printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n"}NR%4==3{print}NR%4==0{printf "##############################\\n"} '  > ''' +
                                 fake_fq2_from_fq1 + "\n")
                fake_fq1_from_fq2 = prefix + "_clean_E1.fq"
                file_shell.write('cat ' + read2_ori + ' | ' + " \\\n" +
                                 ''' awk 'NR%4==1{print}NR%4==2{printf "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n"}NR%4==3{print}NR%4==0{printf "##############################\\n"} '  > ''' +
                                 fake_fq1_from_fq2 + "\n")
                file_shell.write("cat " + read1_rename + ' >> ' + final_out1 + "\n" +
                                 "cat " + fake_fq2_from_fq1 + ' >> ' + final_out2 + "\n" +
                                 "cat " + fake_fq1_from_fq2 + ' >> ' + final_out1 + "\n" +
                                 "cat " + read2_rename + ' >> ' + final_out2 + "\n\n")
    file_shell.close()

#######################################################
#  Main function
#######################################################

if __name__ == '__main__':
    main()
