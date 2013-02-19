#!/usr/bin/env python
# Andrew Borgman
# Takes a folder containing subfolders with tophat alignments created by 
# align_and_quant.py, counts hits for differential expression analysis per individual
# storing information in specified directory

import sys,os

if len(sys.argv > 3):
	print "Usage ./count_hits.py data_folder output_dir"

else:
	# Getting locations of transcript files for each sample
	data_dir = sys.argv[1]
	out_dir = sys.argv[2]

	samples = os.listdir(data_dir)
	try:
		# Needed if this was used to unzip reads into read folders
		samples.remove("unzip_reads.sh")
	except:
		pass

	# Loop through all samples in sub folders
	for s in samples:
		transcript_file = "%s/%s/cufflinks_aligns/transcripts.gtf" % (data_dir, s)
		bam_file = "%s/%s/tophat_aligns/accepted_hits.bam" % (data_dir, s)

		# Outptutting read counts per gene per individual
		ofs = "%s/%s.tsv" % (out_dir, s)
		cmd = "samtools view %s | htseq-count - /share/hg19/annotation/hg19_coding_sequence_only.gtf > %s -s no" % (bam_file, ofs)
		os.system(cmd)
