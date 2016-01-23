# Minnesota Supercomputing Institute
# Tim Nguyen
# 2015

from Bio import SeqIO
import main as master
import time

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

		rotate_handle = open("rotate.txt", "w")	

		for i in range(0, len(sequence), 500):
			wrap_sequence = sequence[i:] + sequence[:i]

			# Two possible methods: wavelet transforms or multiscaling
			if use_wavelet:
				wavelet = master.WaveletTransform(options)
				details_filename = wavelet.dwt(wrap_sequence)
				wavelet_data = wavelet.get_data_from_R(details_filename)
				origin_scaled = wavelet.get_origin(wavelet_data)
				origin = int(origin_scaled * len(sequence))
			else:
				multiscale = master.Multiscaling(options)
				origin = multiscale.multiscaling(wrap_sequence, in_fasta_handle, record_id)

			origin += i
			while origin >= len(sequence):
				origin -= len(sequence)
			distance = min(origin, len(sequence) - origin)
			print "i:", i, ", origin:", origin, ", distance:", distance
			rotate_handle.write(str(i) + "\t" + str(distance) + "\n")

		rotate_handle.close()
		#rearrange_fasta(origin, record_id, sequence, options)

	in_fasta_handle.close()

	print "Runtime:", time.time() - start_time, "seconds"

###############################################################################

if __name__ == '__main__':
	main()