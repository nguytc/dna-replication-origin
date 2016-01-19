# Minnesota Supercomputing Institute
# Tim Nguyen
# 2015

from Bio import SeqIO
from argparse import ArgumentParser
import time
import sys

###############################################################################
# MULTISCALING
###############################################################################
class Multiscaling:
	""" This class contains all the methods required for the multiscaling
		method of obtaining the DNA replication origin. """

	def __init__(self, options):
		self.options = options

	# MULTISCALING ============================================================
	def multiscaling(self, sequence, record_id):
		""" Use the multiscaling method to predict origin. """
		
		subsequence = sequence # constantly shrinking subset of sequence

		seq_offset = 0
		width = self.options.size

		while width >= 128:
			# width cannot be greater than the sequence length
			if width > len(subsequence):
				width /= 2
				continue

			# Determine gap between windows of base pairs given in terms of
		   	# number of base pairs. Will be 1/100 of the window width.
			window_gap = width / 100
			if window_gap < 1:
				window_gap = 1

			skews = self.get_skews(subsequence, width, window_gap)
			slopes = self.get_slopes(skews, width/window_gap)

			if self.options.debug:
				self.add_tab_entries(skews, slopes, width, record_id, seq_offset, sequence)

			# find origin of replication from slopes for current subsequence.
			current_origin = max(slopes, key = slopes.get)

			seq_offset = self.get_seq_offset(current_origin, width, seq_offset, sequence)
			subsequence = self.get_subsequence(subsequence, current_origin, width)

			width /= 2 # halve the window size every iteration
		
		# Determine final replication origin in the ORIGINAL sequence.
		final_origin = current_origin + seq_offset
		if final_origin >= len(sequence):
			final_origin -= len(sequence)

		return final_origin

	# GET SKEWS ============================================================
	def get_skews(self, sequence, width, window_gap, kind = 'GC'):
		"""	Gets skews for each window of base pairs of sequence.
		   	Returns a list of tuples with positions as x and skew values as y.
			NOTE: Indices MUST start at 0 to cooperate with get_slopes function. """
		
		skews = []
		start_index = len(sequence) - width / 2 # Start with center at 0
		at_end = True
		base_frequencies = {}

		while start_index < len(sequence) and at_end or \
			  start_index < len(sequence) - width / 2 and not at_end:
			focus_window = self.get_focus_window(sequence, start_index, width)

			base_frequencies = self.get_base_frequencies(sequence, focus_window, 
								start_index, window_gap, width, base_frequencies)

			# Determine position in sequence that is center of window.
			window_center = start_index + width / 2
			if window_center >= len(sequence):
				window_center -= len(sequence)

			# can specify which skew is desired
			skew = self.calculate_skew(base_frequencies, kind)
			skews.append((window_center, skew))

			start_index += window_gap
			if start_index >= len(sequence):
				start_index -= len(sequence)
				at_end = False

		return skews

	# GET FOCUS WINDOW ======================================================
	def get_focus_window(self, sequence, start_index, width):
		"""	Get the window of base pairs surrounding the previous level's origin
			and return it as a string. """

		end_index = start_index + width

		if end_index < len(sequence):
			focus_window = sequence[start_index:end_index]
		else: # Wrap around to front of fasta file due to circular DNA
			front_length = end_index - len(sequence)
			focus_window = sequence[start_index:len(sequence)] + sequence[:front_length]

		return focus_window

	# GET BASE FREQUENCIES ====================================================
	def get_base_frequencies(self, sequence, focus_window, start_index, window_gap, 
							width, base_frequencies):
		"""	Get the frequencies of the bases in a sequence. """
		
		if start_index == len(sequence) - width / 2:
			# First window requires counting base pairs in entire window
			base_frequencies = self.count_base_frequencies(focus_window)
		else:
			previous_bases = self.get_previous_bases(sequence, start_index, window_gap)
			next_bases = self.get_next_bases(focus_window, width, window_gap)
			base_frequencies = self.adjust_base_frequencies(base_frequencies, previous_bases, next_bases)
		return base_frequencies

	# COUNT FREQUENCY OF BASES ================================================
	def count_base_frequencies(self, sequence):
		"""	Counts the frequency of G, C, A, and T bases in a given sequence.
		   	Returns a dict with bases as keys and frequencies as values. """
		
		base_frequencies = {}
		G = 0.0
		C = 0.0
		A = 0.0
		T = 0.0

		for base in sequence:
			if base == 'G':
				G += 1
			elif base == 'C':
				C += 1
			elif base == 'A':
				A += 1
			elif base == 'T':
				T += 1

		base_frequencies['G'] = G
		base_frequencies['C'] = C
		base_frequencies['A'] = A
		base_frequencies['T'] = T
		return base_frequencies

	# GET PREVIOUS BASES ======================================================
	def get_previous_bases(self, sequence, start_index, window_gap):
		"""	Gets the bases from the previous window that aren't in focus window
		   e.g.	|  		previous window 			|		
								 |			focus window		 |
				| previous_bases |	window overlap	| next_bases | """

		# previous_bases saves the list of bases in the gap between two windows
		if start_index - window_gap >= 0:
			previous_bases = sequence[start_index - window_gap:start_index]
		else:
			end_length = window_gap - start_index
			previous_bases = sequence[len(sequence) - end_length:len(sequence)] + \
							 sequence[:start_index]
			
		return previous_bases

	# GET NEXT BASES ==========================================================
	def get_next_bases(self, focus_window, width, window_gap):
		"""	Gets the bases from the focus window that weren't in previous window
		   e.g.	|  		previous window 			|		
								 |			focus window		 |
				| previous_bases |	window overlap	| next_bases |"""

		# next_bases is list of bases that weren't considered in last window
		next_bases = focus_window[width - window_gap:width]
		return next_bases

	# ADJUST BASE FREQUENCY COUNT =============================================
	def adjust_base_frequencies(self, base_frequencies, previous_bases, next_bases):
		"""	Adjusts the frequency counts of G, C, A, and T bases for new bases.
		   	Returns a dict with bases as keys and frequencies as values. """

		for previous_base in previous_bases:
			if previous_base == 'G' or previous_base == 'C' or \
			   previous_base == 'A' or previous_base == 'T':
				base_frequencies[previous_base] -= 1
		
		for next_base in next_bases:
			if next_base == 'G' or next_base == 'C' or \
			   next_base == 'A' or next_base == 'T':
				base_frequencies[next_base] += 1

		return base_frequencies

	# CALCULATE SKEW =======================================================
	def calculate_skew(self, base_frequencies, kind = 'GC'):
		"""	Calculates either the GC, AT, or GT-AC skew given the frequencies 
			of the bases. """

		G = base_frequencies['G']
		C = base_frequencies['C']
		A = base_frequencies['A']
		T = base_frequencies['T']

		if kind == 'GC' and (G + C) != 0: 
			skew = (G - C) / (G + C)
		elif kind == 'AT' and (A + T) != 0:
			skew = (A - T) / (A + T)
		elif kind == 'GTAC' and (A*G + A*C + T*G + T*C) != 0:
			skew = (2*A*G - 2*T*C) / (A*G + A*C + T*G + T*C)
		else:
			print "ERROR: calculate_skew"
			skew = 0

		return skew

	# GET SLOPES ==============================================================
	def get_slopes(self, skews, skews_window_size):
		"""	Gets the skew slopes of positions given by linear regression.
		   	Returns a dict with positions as key and slopes as value. """

		slopes = {}
		start_index = 0
		# There will be 2 entries with an index of skews[-1][0] due to 0 at start
		last_index_taken = False

		while start_index < len(skews):
			skews_tuples = self.get_skews_tuples(skews, start_index, skews_window_size)
			slope = self.calculate_slope(skews_tuples)

			# skews_window_center is the position given by the middle tuple
			skews_window_center = skews_tuples[(len(skews_tuples) - 1) / 2][0]

			if skews_window_center > skews[-1][0] or \
			   skews_window_center == skews[-1][0] and last_index_taken:
				skews_window_center -= skews[-1][0]
			elif skews_window_center == skews[-1][0] and not last_index_taken:
				last_index_taken = True
				
			slopes[skews_window_center] = slope

			start_index += 1
			
		return slopes

	# GET POSITION, SKEW VALUE TUPLES =========================================
	def get_skews_tuples(self, skews, start_index, skews_window_size):
		"""	Get the window of position and skew value tuples.
		   	Returns the tuples in a list. """

		end_index = start_index + skews_window_size

		if end_index < len(skews):
			skews_tuples = skews[start_index:end_index]
		else: 
			# Wrap around to front of skews list because of circular DNA
			# The indices of the front pairs must be adjusted for a correct slope
			front_length = end_index - len(skews)
			adjust_amount = skews[-1][0] # Add index of last element to each front
			front_adjusted = [(pair[0] + adjust_amount, pair[1]) for pair in skews[:front_length]]
			skews_tuples = skews[start_index:len(skews)] + front_adjusted

		return skews_tuples

	# CALCULATE LINEAR REGRESSION SLOPE =======================================
	def calculate_slope(self, xy_tuples):
		""" Calculates the slope of the linear regression line of the data set.
		    xy_tuples is a set of x and y coordinates.
		    m = (n*sum(xy)-sum(x)*sum(y)) / (N*sum(x^2) - sum(x)^2) """

		n = len(xy_tuples)
		x_sum = 0.0
		y_sum = 0.0
		xy_sum = 0.0
		x_squared_sum = 0.0

		for xy_tuple in xy_tuples:
			x = xy_tuple[0]
			y = xy_tuple[1]
			x_sum += x
			y_sum += y
			xy_sum += x*y
			x_squared_sum += x**2

		if (n*x_squared_sum - x_sum**2) != 0:
			slope = (n*xy_sum - x_sum*y_sum) / (n*x_squared_sum - x_sum**2)
		else:
			print "ERROR: calculate_slope divide by 0"
			slope = 0

		return slope

	# GET SEQUENCE OFFSET =====================================================
	def get_seq_offset(self, origin, width, seq_offset, sequence):
		""" Keeps track of which index the current subsequence starts at in
			relation to the original sequence. Returns the start index. """

		seq_offset += origin - width / 2
		if seq_offset < 0:
			seq_offset += len(sequence) # len(sequence) + (-x)
		elif seq_offset >= len(sequence):
			seq_offset -= len(sequence)
		return seq_offset

	# GET SUBSEQUENCE =========================================================
	def get_subsequence(self, sequence, origin, width):
		"""	Given an origin and width, obtain window surrounding origin. """

		begin_position = origin - width / 2
		end_position = begin_position + width

		if begin_position >= 0 and end_position < len(sequence):
			subsequence = sequence[begin_position:end_position]
		elif begin_position < 0 and end_position < len(sequence):
			back_begin_position = len(sequence) + begin_position # len(sequence) + (-x)
			subsequence = sequence[back_begin_position:len(sequence)] + sequence[:end_position]
		elif begin_position >= 0 and end_position >= len(sequence):
			front_end_position = end_position - len(sequence)
			subsequence = sequence[begin_position:len(sequence)] + sequence[:front_end_position]
		else: # This should never happen
			print "ERROR: get_subsequence window width larger than sequence"
			subsequence = ""

		return subsequence

	# ADD ENTRIES TO TAB-DELIMITED OUTPUT FILE ================================
	def add_tab_entries(self, skews, slopes, width, record_id, seq_offset, sequence):
		"""	Write GC skew values to the tab-delimited output file.
		   	Each row contains a position, window, value, slope, and record ID. """
		
		tab_out_handle = open(self.options.out_filename + '_output.tab', 'w')

		row_entry_list = ['', '', '', '', '']
		for skew in skews:
			position = skew[0]
			skew_value = skew[1]

			row_entry_list[0] = position + seq_offset
			if row_entry_list[0] >= len(sequence):
				row_entry_list[0] -= len(sequence)

			row_entry_list[1] = width
			row_entry_list[2] = skew_value
			row_entry_list[3] = slopes[position]
			row_entry_list[4] = '\"' + record_id + '\"'
			tab_out_handle.write('\t'.join(str(x) for x in row_entry_list) + '\n')

		tab_out_handle.close()

###############################################################################
# WAVELET TRANSFORM
###############################################################################
class WaveletTransform:
	""" This class contains all the methods required for the wavelet transform
		method of obtaining the DNA replication origin. 
		More information: http://www.pybytes.com/pywavelets/ """

	def __init__(self, options):
		# modules required for wavelet transform method but not for multiscaling
		from rpy2 import robjects
		import pywt

		self.data = {}
		self.rank_lengths = {}
		self.options = options

	class RankSegment:
		def __init__(self, min_pos, max_pos, coeff):
			self.children = []
			self.min_pos = min_pos
			self.max_pos = max_pos
			self.coeff = coeff

	# WAVELET TRANSFORM =======================================================
	def wavelet_transform(self,sequence):
		""" Use the wavelet transform method to predict the origin. """

		details_filename = self.dwt(sequence)
		wavelet_data = self.get_data_from_R(details_filename)
		origin_scaled = self.get_origin(wavelet_data)
		final_origin = int(origin_scaled * len(sequence))

		return final_origin

	# DISCRETE WAVELET TRANSFORM ==============================================	
	def dwt(self, sequence):
		""" Use multilevel decomposition Haar discrete wavelet transform method 
			with periodic-padding. 
			Outputs a file with the details coefficients. """

		details_filename = self.options.out_filename + '_details.txt'
		details_handle = open(details_filename, 'w')

		# encode sequence of ATCG characters to be integers from -1 to 1
		sequence_encoding = self.encode_sequence(sequence)

		w = pywt.Wavelet('haar')
		max_level = pywt.dwt_max_level(len(sequence_encoding), w)

		coeffs = pywt.wavedec(sequence_encoding, 'haar', 'ppd', level=max_level)
		cA_n = coeffs[0] # approximation coefficients not needed
		cD_arrays = coeffs[1:]

		for cD in cD_arrays:
			details_handle.write('\t'.join(str(d) for d in cD) + '\n')

		details_handle.close()
		return details_filename

	# ENCODE SEQUENCE =========================================================
	def encode_sequence(self, sequence):
		""" Given a string of ACTG characters representing base pairs, encode them 
			to be array of integers for wavelet transform.
			G = -1, A = 0, T = 0, C = 1. Return encoded array. """

		sequence_encoding = []
		for base in sequence:
			if base == 'G':
				sequence_encoding.append(-1)
			elif base == 'C':
				sequence_encoding.append(1)
			else: # base is A, T, or other
				sequence_encoding.append(0)
		return sequence_encoding

	# GET WAVELET DATA FROM R =================================================
	def get_data_from_R(self, details):
		""" Obtains quantitative data for the wavelet transform from R. """

		# body of the R script
		robjects.r("""
		plot <- function(details) {
			library(ggplot2)
			test = scan(details, what='character', sep="\\n")
			test2 = lapply(test, function(x){return(unlist(strsplit(x,"\\t")))})
			wavelet_data = data.frame()
			for(i in 1:length(test2)){
	 			SiteCount = length(test2[[i]])
	 			SiteIntervals = seq(0,1,by=1/SiteCount)
	 			wavelet_data= rbind(wavelet_data, data.frame(Rank=i, PosMin=SiteIntervals[-1*length(SiteIntervals)], PosMax=SiteIntervals[-1], Coeffs=as.numeric(test2[[i]])/sd(as.numeric(test2[[i]]))))
			}
			return(wavelet_data)
		}
		""")

		plot_f = robjects.r['plot']	
		wavelet_data = zip(*plot_f(details))
		return wavelet_data

	# GET REPLICATION ORIGIN USING MULTIPLE PATHS =============================
	def get_origin(self, wavelet_data):
		""" Finds the origin for DNA replication given the data from wavelet. 
			Returns the position of the origin scaled within the sequence as a 
			number from 0.0 to 1.0 """

		self.put_data_in_dict(wavelet_data)
		self.find_rank_lengths()
		self.determine_children()
		(coeff_sum, origin_segment) = self.find_origin_path((0, 0), 0)

		origin_scaled = (self.data[origin_segment].min_pos + \
						 self.data[origin_segment].max_pos) / 2
		return origin_scaled

	# CONVERT WAVELET DATA ARRAY TO DICTIONARY ================================
	def put_data_in_dict(self, wavelet_data):
		""" Given a plot array containing the wavelet transform data, store the
			information into an easily accessible data structure.
			Stores data in a dictionary with keys being a tuple (rank, segment) 
			and the corresponding values being a list containing 3 elements: 
			min pos, max pos, and details coefficient. 
			i.e. (rank, segment), [min pos, max pos, details coefficient] """
  
		# e.g., | segment 1 | segment 2	| segment 3	| segment 4	|
		rank = 1 # level in the wavelet transform
		segment = 1 # position in the current level  

		for row in wavelet_data:
			if row[0] != rank:
				rank = row[0]
				segment = 1

			rank_segment = self.RankSegment(row[1], row[2], row[3])
			self.data[(rank, segment)] = rank_segment
			segment += 1

	# FIND LENGTHS OF RANKS IN TERMS OF SEGMENTS ==============================
	def find_rank_lengths(self):
		""" Get the keys of the data dictionary and organize them by ranks.
			Creates a dictionary with keys as ranks and values as the number of
			segments in the data dictionary in that rank. 
			e.g., (1, 5) means there are 5 segments in rank 1 """

		for segment in self.data:
			rank = segment[0]

			if rank not in self.rank_lengths:
				self.rank_lengths[rank] = 1
			else:
				self.rank_lengths[rank] += 1

	# DETERMINE CHILDREN ======================================================
	def determine_children(self):
		""" Given a dictionary containing the wavelet transform data, determine 
			all the children of each rank segment in the wavelet transform.
			Stores the children as a list of tuples that may be reached from
			the key position. 
			e.g., (rank + 1, 0), (rank + 1, 1), (rank + 1, 2) """

		for segment in self.data:
			rank = segment[0]
			if rank == len(self.rank_lengths):
				continue

			children = self.data[segment].children
			
			start = 1
			end = self.rank_lengths[rank + 1]
			focus_min_position = self.data[segment].min_pos
			focus_max_position = self.data[segment].max_pos
			
			# find where possible children could start with binary search
			while start <= end:
				mid = (start + end) / 2
				child_min_position = self.data[(rank + 1, mid)].min_pos
				child_max_position = self.data[(rank + 1, mid)].max_pos

				if child_min_position <= focus_max_position and \
				   child_max_position >= focus_min_position:
					break
				elif child_min_position > focus_max_position:
					end = mid - 1
				elif child_max_position < focus_min_position:
					start = mid + 1

			# back the mid up to include all possible children before the mid
			while mid > 1 and child_max_position > focus_min_position:
				mid -= 1
				child_max_position = self.data[(rank + 1, mid)].max_pos

			child_min_position = self.data[(rank + 1, mid)].min_pos
			while child_min_position < focus_max_position:
				if child_min_position < focus_max_position and \
				   child_max_position > focus_min_position:
					children.append((rank + 1, mid))

				mid += 1
				if mid > self.rank_lengths[rank + 1]:
					break
				child_min_position = self.data[(rank + 1, mid)].min_pos
				child_max_position = self.data[(rank + 1, mid)].max_pos

	# FIND PATH WITH MAX COEFFICIENTS SUM =====================================
	def find_origin_path(self, segment, coeff_sum):
		""" Obtain the origin from the wavelet data using DFS with recursion. 
			Sum up the coefficients of each path to determine the path with the 
			largest sum of coefficients. """

		# max of the summed up coefficients
		max_coeff_sum = None
		# rank and segment of the curent origin (rank, segment)
		origin = None

		# start of the tree. beginning of the recursion.
		if segment == (0, 0):
			(sum_1, origin_1) = self.find_origin_path((1, 1), self.data[(1, 1)].coeff)
			(sum_2, origin_2) = self.find_origin_path((1, 2), self.data[(1, 2)].coeff)
			if sum_2 > sum_1:
				max_coeff_sum = sum_2
				origin = origin_2
			else:
				max_coeff_sum = sum_2
				origin = origin_1

		# end of the recursion at the rank second to the lowest
		elif segment[0] == len(self.rank_lengths) - 1:
			# determine which child at the bottom has largest coefficient
			for child in self.data[segment].children:
				now_coeff_sum = coeff_sum + self.data[child].coeff
				if max_coeff_sum == None or now_coeff_sum > max_coeff_sum:
					max_coeff_sum = now_coeff_sum
					origin = child

		else:
			for child in self.data[segment].children:
				now_coeff_sum = coeff_sum + self.data[child].coeff
				(focus_sum, focus_origin) = self.find_origin_path(child, now_coeff_sum)
				if max_coeff_sum == None or focus_sum > max_coeff_sum:
					max_coeff_sum = focus_sum
					origin = focus_origin

		return (max_coeff_sum, origin)

	# STATIONARY WAVELET TRANSFORM ============================================
	def _swt(self, sequence):
		""" Use stationary wavelet transform method. """

		details_filename = self.options.out_filename + '_details.txt'
		details_handle = open(details_filename, 'w')

		padded_sequence = self._pad_sequence(sequence)
		sequence_encoding = self.encode_sequence(padded_sequence)
		level = pywt.swt_max_level(len(sequence_encoding))
		coeffs = pywt.swt(sequence_encoding, 'haar', level)

		for (cA, cD) in coeffs:
			approx_handle.write('\t'.join(str(a) for a in cA) + '\n')
			details_handle.write('\t'.join(str(d) for d in cD) + '\n')

		details_handle.close()
		return details_filename

	# PAD SEQUENCE ============================================================
	def _pad_sequence(self, sequence):
		""" Pad given sequence to have a length of a power of two.
			Pad with the beginning of the sequence due to circular DNA. """

		current_power_of_2 = 1
		while (len(sequence) > current_power_of_2):
			current_power_of_2 *= 2

		if len(sequence) == current_power_of_2:
			padded_sequence = sequence
		else:
			number_of_extra = current_power_of_2 - len(sequence)
			padding = sequence[:number_of_extra]
			padded_sequence = sequence + padding
			
		return padded_sequence

	# GET REPLICATION ORIGIN USING SINGLE PATH ================================
	def _get_origin_single_path(self, wavelet_data):
		""" Travel down wavelet data through the highest coefficient of each
			level. 
			NOTE: This function is not used. get_origin is used instead """

		# master positions cannot be overstepped in the sublevels
		master_min_position = 0.0
		master_max_position = 1.0

		# regular positions are from the level currently being considered
		min_position = master_min_position
		max_position = master_max_position
	 
		focus_rank = 1 # rank is the level in the wavelet transform
		level_high_coeff = -sys.maxint - 1

		# row[0] = rank, row[1] = min pos, row[2] = max pos, row[3] = coeff
		for row in wavelet_data:
			# the level has changed
			if row[0] != focus_rank:
				focus_rank = row[0]
				if min_position > master_min_position:
					master_min_position = min_position
				if max_position < master_max_position:
					master_max_position = max_position
				level_high_coeff = -sys.maxint - 1
				min_position = master_min_position
				max_position = master_max_position

			if row[1] > master_max_position or row[2] < master_min_position:
				continue
			elif row[3] > level_high_coeff:
				min_position = row[1]
				max_position = row[2]
				level_high_coeff = row[3]
		
		if min_position > master_min_position:
			master_min_position = min_position
		if max_position < master_max_position:
			master_max_position = max_position		
		return (master_min_position + master_max_position) / 2.0

###############################################################################
# Z CURVE
###############################################################################
class ZCurve:
	""" This class contains all the methods required for the z curve
		method of obtaining the DNA replication origin. 
		NOTE: Not fully implemented. The data from the z curve method is
		obtained but the data is not analyzed to get the origin. """

	def __init__(self, options):
		self.options = options

	# ADD ENTRIES TO TAB-DELIMITED OUTPUT FILE ================================
	def _add_tab_entries(self, n, x, y, z, result, tab_out_handle, xy_handle):
		"""	Write values to the tab-delimited output file.
		   	Each row contains a position, x, y, z, and result. """

		row_entry_list = ['', '', '', '', '']

		row_entry_list[0] = n
		row_entry_list[1] = x
		row_entry_list[2] = y
		row_entry_list[3] = z
		row_entry_list[4] = result
		tab_out_handle.write('\t'.join(str(x) for x in row_entry_list) + '\n')

		xplusy = x + y
		xy_handle.write(str(n) + '\t' + str(xplusy) + '\n')

	# Z CURVE =================================================================
	def _zcurve(self, sequence):
		""" Use z curve method and calculate x, y, and z values of positions. 
			Outputs results in tab-delimited file and new fasta file. """

		tab_out_handle = open('zcurve_output.tab', 'w')
		xy_handle = open('zcurve_xplusy.tab', 'w')

		G = 0.0
		C = 0.0
		A = 0.0
		T = 0.0
		n = 1.0

		for base in sequence:
			if base == 'G':
				G += 1
			elif base == 'C':
				C += 1
			elif base == 'A':
				A += 1
			elif base == 'T':
				T += 1

			x = (A + G) - (C + T)
			y = (A + C) - (G + T)
			z = (A + T) - (G + C)

			if (n - z) != 0:
				result = (y - x) / (n - z)
			else:
				result = 0

			self._add_tab_entries(n, x, y, z, result, tab_out_handle, xy_handle)

			n += 1

		tab_out_handle.close()
		xy_handle.close()

# PARSE COMMAND LINE ARGUMENTS ------------------------------------------------
def parse_arguments():
	"""	Parses the given command line options and arguments
	   	Returns a dictionary with the given arguments as the values. """	
	
	use = "python origin.py -f <file.fasta> [options]"
	parser = ArgumentParser(usage=use)

	parser.add_argument('-f', '--fasta', dest='in_fa_filename', 
		help = "Required input FASTA (.fasta/.fa) file")
	parser.add_argument('-o', '--out', dest='out_filename', default='adjusted',
		help = "Basename of rearranged FASTA output file [adjusted]")
	parser.add_argument('-s', '--size', dest='size', type=int, default=524288, 
		help = "Starting window size. Decreases by 2 every iteration [524288]")
	parser.add_argument('-d', '--debug', dest='debug', action='store_true', 
		default=False, help = "Prints out all intermediate files obtained [False]")

	# Specify which method is desired
	parser.add_argument('-m', '--multiscale', dest='use_multiscale', action='store_true', 
		default=True, help = "Use the multiscaling method (Recommended) [True]")
	parser.add_argument('-w', '--wavelet', dest='use_wavelet', action='store_true', 
		default=False, help = "Use the wavelet transform method instead [False]")

	options = parser.parse_args()	

	# FASTA file REQUIRED. Display error and exit if not given in command line
	if not options.in_fa_filename:
		parser.print_help()
		print "\nERROR: Fasta file not given"
		sys.exit()

	return options

# REARRANGE FASTA FILE --------------------------------------------------------
def rearrange_fasta(origin, record_id, sequence, options):
	"""	Rearrange sequence to begin at DNA replication origin and loops around.
	   	Writes the rearranged sequence to a fasta file. """

	out_fasta_handle = open(options.out_filename + '_output.fa', 'w')
	rearranged_sequence_newlines = [] # rearranged sequence with newlines added
	
	out_fasta_handle.write('>' + record_id + '\n')
	rearranged_sequence = sequence[origin:] + sequence[:origin]

	# Insert a newline every 70 base pairs
	for i in xrange(0, len(rearranged_sequence), 70):
		rearranged_sequence_newlines.append(rearranged_sequence[i: i + 70])
	out_fasta_handle.write('\n'.join(rearranged_sequence_newlines) + '\n')

	out_fasta_handle.close()

# MAIN ------------------------------------------------------------------------
def main():
 	start_time = time.time()

	options = parse_arguments()
	in_fasta_handle = open(options.in_fa_filename, 'rU')
	
	# script not fully tested with multiple records
	for record in SeqIO.parse(in_fasta_handle, 'fasta'):
		record_id = record.id
		sequence = str(record.seq)
		print "Sequence Length:", len(sequence)

		# Two possible methods: wavelet transforms or multiscaling
		if options.use_wavelet:
			wavelet = WaveletTransform(options)
			origin = wavelet.wavelet_transform(sequence)
		else:
			multiscale = Multiscaling(options)
			origin = multiscale.multiscaling(sequence, record_id)

		print "DNA Replication Origin:", origin
		rearrange_fasta(origin, record_id, sequence, options)

		print "Rearranged FASTA file:", options.out_filename + "_output.fa"

	in_fasta_handle.close()

	print "Runtime:", time.time() - start_time, "seconds"

###############################################################################

if __name__ == '__main__':
	main()
