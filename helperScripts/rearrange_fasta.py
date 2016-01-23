# Minnesota Supercomputing Institute
# Tim Nguyen
# 2015

from Bio import SeqIO
import time
import sys

# script used to rearrange fasta file by however many bases each time
# this is a scrappy script

# REARRANGE FASTA FILE --------------------------------------------------------
def rearrange_fasta(origin, record_id, sequence):
	"""	Rearrange sequence to begin at DNA replication origin and loops around.
	   	Writes the rearranged sequence to a fasta file. """

	out_fasta_handle = open(str(origin) + '_output.fa', 'w')
	rearranged_sequence_newlines = [] # rearranged sequence with newlines added
	
	out_fasta_handle.write('>' + record_id + '\n')
	rearranged_sequence = sequence[origin:] + sequence[:origin]

	# Insert a newline every 70 base pairs
	for i in xrange(0, len(rearranged_sequence), 70):
		rearranged_sequence_newlines.append(rearranged_sequence[i: i + 70])
	out_fasta_handle.write('\n'.join(rearranged_sequence_newlines) + '\n')

	out_fasta_handle.close()

# MAIN ------------------------------------------------------------------------
def main(argv):
 	start_time = time.time()

	in_fasta_handle = open(argv[1], 'rU')
	
	print argv[1]

	# script not fully tested with multiple records
	for record in SeqIO.parse(in_fasta_handle, 'fasta'):
		record_id = record.id
		sequence = str(record.seq)

		for i in range(0, len(sequence), 200000):

			print "DNA Replication Origin:", i
			rearrange_fasta(i, record_id, sequence)

			print "Rearranged FASTA file:" + str(i) + "_output.fa"

	in_fasta_handle.close()

	print "Runtime:", time.time() - start_time, "seconds"

###############################################################################

if __name__ == '__main__':
	main(sys.argv)