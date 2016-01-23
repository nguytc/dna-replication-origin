# Minnesota Supercomputing Institute
# Tim Nguyen
# 2015

from Bio import SeqIO
from random import randrange
import main as master
import time

# MUTATE SEQUENCE -------------------------------------------------------------
def mutate_sequence(rate, sequence):
	""" Randomly mutates a given sequence with a specified mutation rate.
		rate -> int representing percent. eg. 20 is 20 pct mutation rate. 
		Returns a mutated sequence as a string. Used for testing purposes. """

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
		
	return "".join(mutated_sequence)

# MAIN ------------------------------------------------------------------------
def main():
 	start_time = time.time()

	options = master.parse_arguments()
	use_wavelet = options.wavelet
	in_fasta_handle = open(options.in_fasta_filename, "rU")
	
	# script not fully tested with multiple records
	for record in SeqIO.parse(in_fasta_handle, "fasta"):
		record_id = record.id
		sequence = str(record.seq)

		mutate_handle = open("mutate.txt", "w")	

		for rate in range(1, 101):
			for i in range(100):
				mutated_sequence = mutate_sequence(rate, sequence)

				# Two possible methods: wavelet transforms or multiscaling
				if use_wavelet:
					wavelet = master.WaveletTransform(options)
					details_filename = wavelet.dwt(mutated_sequence)
					wavelet_data = wavelet.get_data_from_R(details_filename)
					origin_scaled = wavelet.get_origin(wavelet_data)
					origin = int(origin_scaled * len(sequence))
				else:
					multiscale = master.Multiscaling(options)
					origin = multiscale.multiscaling(mutated_sequence, in_fasta_handle, record_id)

				distance = min(origin, len(sequence) - origin)
				print "rate:", rate, ", i:", i, ", distance:", distance
				mutate_handle.write(str(rate) + "\t" + str(distance) + "\n")

		mutate_handle.close()
		#rearrange_fasta(origin, record_id, sequence, options)

	in_fasta_handle.close()

	print "Runtime:", time.time() - start_time, "seconds"

###############################################################################

if __name__ == '__main__':
	main()