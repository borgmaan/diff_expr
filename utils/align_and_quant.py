#!/usr/bin/env python
# Andrew Borgman
# 1/24/2013
# Takes a file listing sample names and folder locations
# Aligns reads against ref using tophat and quatifies 
# trascripts using cufflinks

import sys,os
import csv

def analyze_run(run_folder="",home_directory=""):
    os.chdir(run_folder)
    os.system("rm -r tophat_all")
    os.system("rm -r cufflinks_all")
    os.system("mkdir tophat_all")
    os.system("mkdir cufflinks_all")
    tophat_output_loc = os.getcwd() + "/tophat_all/"
    cufflinks_output_loc = os.getcwd() + "/cufflinks_all/"
    tophat_string = " ".join(["tophat2","-o",tophat_output_loc,"-p 10", "-r 250", "--b2-very-sensitive","/share/hg19/bowtie_indexes/hg19","reads_1.fastq","reads_2.fastq"])
    cuff_string = " ".join(["cufflinks","-o",cufflinks_output_loc,"-p 10", "-G /share/hg19/hg19.gtf",tophat_output_loc + "accepted_hits.bam"])
    os.system(tophat_string)
    os.system(cuff_string)
    os.chdir(home_directory)
    return None


if __name__ == "__main__":
    home_dir = os.getcwd()
    zip_dir = "/home/andrew/bernie_analysis/data_run_4/Data/"
    with open("run_4_map") as mapper:
        for line in mapper:
            spl = [x.replace("\n","").strip() for x in line.split()]
            print spl

            sub_folder = spl[0]
            sample_name = spl[1]
            read_file_1 = spl[2].split(",")[0]
            read_file_2 = spl[2].split(",")[1]
            sample_dir = "./data_run_4/%s/" % sub_folder
            print "Making:",sub_folder
            os.system("mkdir " + sample_dir) 
            zip_string_1 = "gunzip -c %s%s > %sreads_1.fastq" %(zip_dir, read_file_1, sample_dir)
            zip_string_2 = "gunzip -c %s%s > %sreads_2.fastq" %(zip_dir, read_file_2, sample_dir)
            print "Unzipping files"
            #os.system(zip_string_1) 
            #os.system(zip_string_2) 
            print "Analyzing"
            analyze_run(sample_dir,home_dir) 
