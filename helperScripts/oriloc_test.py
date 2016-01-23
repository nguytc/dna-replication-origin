# Minnesota Supercomputing Institute
# Tim Nguyen
# 2015

from rpy2 import robjects
from Bio import SeqIO
from random import randrange
import time
import argparse
import sys

def main():
	start_time = time.time()

	options = parse_arguments()
	use_rotate = options.rotate
	use_mutate = options.mutate

	robjects.r("""
	plot <- function(fasta, predict) {
		library(seqinr)
		out = oriloc(seq.fasta = fasta, g2.coord = predict)
		row = which.max(out[,5])
		row_data = out[row,]

		start_b = row_data[,2] * 1000
		end_b = row_data[,3] * 1000
		origin = (start_b + end_b) / 2
		return(origin)
	}
	""")

	plot_f = robjects.r["plot"]

	in_fasta_handle = open(options.in_fasta_filename, "rU")

	for record in SeqIO.parse(in_fasta_handle, "fasta"):
		record_id = record.id
		sequence = str(record.seq)

		if use_rotate:
			# Rotations =======================================================
			rotate_handle = open("rotate_oriloc.txt", "w")	
			for offset in range(0, len(sequence), 500):
				rearrange_fasta(offset, record_id, sequence)
				adjust_predict(offset, sequence, options.in_predict_filename)

				origin = int(plot_f("ori_test.fa", "ori_test.predict")[0]) + offset
				while origin >= len(sequence):
					origin -= len(sequence)
				distance = min(origin, len(sequence) - origin)
				print "offset:", offset, ", distance:", distance
				rotate_handle.write(str(offset) + "\t" + str(distance) + "\n")

			rotate_handle.close()

		if use_mutate:
			# Mutations ===========================================================
			mutate_handle = open("mutate_oriloc.txt", "w")
			for rate in range(1, 101):
				for i in range(100):
					mutate_sequence(rate, record_id, sequence)

					origin = int(plot_f("ori_test.fa", options.in_predict_filename)[0])
					distance = min(origin, len(sequence) - origin)
					print "rate:", rate, ", i:", i, ", distance:", distance
					mutate_handle.write(str(rate) + "\t" + str(distance) + "\n")
			mutate_handle.close()
		

	in_fasta_handle.close()

	print "Runtime:", time.time() - start_time, "seconds"

def rearrange_fasta(offset, record_id, sequence):
	out_fasta_handle = open("ori_test.fa", "w")
	rearranged_sequence_newlines = [] # rearranged sequence with newlines added
	
	out_fasta_handle.write(">" + record_id + "\n")
	rearranged_sequence = sequence[offset:] + sequence[:offset]

	# Insert a newline every 70 base pairs
	for i in xrange(0, len(rearranged_sequence), 70):
		rearranged_sequence_newlines.append(rearranged_sequence[i: i+70])
	out_fasta_handle.write('\n'.join(rearranged_sequence_newlines) + "\n")
	out_fasta_handle.close()

def adjust_predict(offset, sequence, predict_filename):
	in_predict_handle = open(predict_filename, "r+")
	out_predict_handle = open("ori_test.predict", "w")

	first_line = True
	for line in in_predict_handle:
		if first_line == True:
			out_predict_handle.write(line)
			first_line = False
			continue

		str_split = line.split( )
		start = int(str_split[1]) + offset
		end = int(str_split[2]) + offset

		while start >= len(sequence):
			start -= len(sequence)
		while end >= len(sequence):
			end -= len(sequence)

		str_split[1] = start
		str_split[2] = end

		out_predict_handle.write(" ".join(str(x) for x in str_split) + "\n")
	
	in_predict_handle.close()
	out_predict_handle.close()

def mutate_sequence(rate, record_id, sequence):
	""" Randomly mutates a given sequence with a specified mutation rate.
		rate -> int representing percent. eg. 20 is 20 pct mutation rate. """

	bases = ['A', 'T', 'C', 'G']

	if rate <= 0:
		return sequence
	elif rate > 100:
		rate = 100

	mutated_sequence = list(sequence)

	for i in range(len(mutated_sequence)):
		mutate = randrange(1, 101)
		if mutate <= rate: # mutate base pair
			while 1:
				mutated_base = bases[randrange(0, 4)]
				if mutated_base != mutated_sequence[i]:
					break
			mutated_sequence[i] = mutated_base
		
	mutated_string = "".join(mutated_sequence)

	mutated_sequence_newlines = []
	out_fasta_handle = open("ori_test.fa", "w")
	out_fasta_handle.write(">" + record_id + "\n")

	# Insert a newline every 70 base pairs
	for i in xrange(0, len(mutated_string), 70):
		mutated_sequence_newlines.append(mutated_string[i: i+70])
	out_fasta_handle.write('\n'.join(mutated_sequence_newlines) + "\n")
	out_fasta_handle.close()

def parse_arguments():
	"""	Parses the given command line options and arguments
	   	Returns a dictionary with the given arguments as the values. """	

	usage = "Usage: python oriloc_test.py -f <file.fasta> [options]"
	parser = argparse.ArgumentParser(usage)

	parser.add_argument("-f", dest="in_fasta_filename", 
		help = "Required FASTA (.fasta/.fa) file")
	parser.add_argument("-p", dest="in_predict_filename", 
		help = "Required predict (.predict) file")
	parser.add_argument("-r", dest="rotate", action="store_true", default=False, 
		help = "Use the rotation testing method [False]")
	parser.add_argument("-m", dest="mutate", action="store_true", default=False, 
		help = "Use the mutation testing method [False]")

	options = parser.parse_args()	

	# FASTA file REQUIRED. Display error and exit if not given in command line
	if not options.in_fasta_filename:
		parser.print_help()
		print "\nERROR: Fasta file not given"
		sys.exit()
	# predict file REQUIRED. Display error and exit if not given in command line
	if not options.in_predict_filename:
		parser.print_help()
		print "\nERROR: Predict file not given"
		sys.exit()

	return options

if __name__ == "__main__":
	main()