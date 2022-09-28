#############################################################################
#Author: Christina Castro
#Date: March 31, 2016
#Affiliation: Centers for Disease Control and Prevention
#			  Division of Viral Diseases
#			  Advanced Molecular Detection
#############################################################################

#USAGE: python U50.py <reference fasta file> <sorted BED file>

#**** ATTENTION: All coordinates are based on the BED file's ZERO-based counting system ****#

import sys
import numpy as np
# np.set_printoptions(threshold="nan")
from Bio import SeqIO


#chromStart - The starting position of the feature in the chromosome or scaffold. 
#			  The first base in a chromosome is numbered 0.
#chromEnd - The ending position of the feature in the chromosome or scaffold. 
#			The chromEnd base is not included in the display of the feature. 
#			*For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, 
#			and span the bases numbered 0-99.



#-----------------------------------------[ Declare Input Variables ]------------------------------------------#

#Input is a single reference fasta file
fasta_file = sys.argv[1]
#Input is a sorted BED file of contigs
BED_file = sys.argv[2]

#------------------------------------------------[ Functions ]-------------------------------------------------#

#####################################################################
#	Loop through FASTA file - Get reference length and create array
#####################################################################
def getfasta(file):
	fasta_sequences = SeqIO.parse(open(file),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
			
		fasta_length = len(sequence)			
		length_array = np.zeros(len(sequence), dtype=np.int8)		#Creates an array of zeros the length of the input reference sequence	
		return fasta_length, length_array

##################################################################################################################
#	Loop through BED file - Assign unique start and stop coordinates and contig length to lists and a dictionary
##################################################################################################################
def getBED(file):
	L = ""
	start_coordinates = []
	stop_coordinates = []
	length_list = []
	contig_dict = {}	
	with open(file, "r")as infile_bed:
		for line in infile_bed:
			length = 0
			L = line.strip()	
			ref, start, stop, read_name, score, strand = L.split("\t")
			stop = int(stop)-1		#Accounts for the BED file format, includes start up to but not including stop
			stop = str(stop)		#Convert back to string
			contig_length = int(stop) - int(start)	
			if (start not in start_coordinates) and (stop not in stop_coordinates):		#Removes any duplicates
				start_coordinates.append(start)
				stop_coordinates.append(stop)
				length_list.append(contig_length)
				contig_dict[str(int(start))+".."+str(int(stop))] = contig_length
		return start_coordinates, stop_coordinates, length_list, contig_dict

####################################################################################
#	Reverse a list
####################################################################################
def reverseList(list):	
	return list.reverse()

####################################################################################
#	Find the overlapping base positions (starts with longest contig to shortest). 
#	Utilizes the array of ones created from fasta file. 
####################################################################################	
def findOverlaps(sorted_start, sorted_stop, length_array):
	Overlapping_Index = {}		
	Unique_Count = {}
	counter = int(0)
	#Checks each contig (starting with the longest) to find overlapping regions
	for i in range(0,len(sorted_start)):
		a = sorted_start[i]
		b = sorted_stop[i]
		if int(b) >= int(len(length_array)):
			b = int(len(length_array) - 1)
		#Loop through the range specified by the start and stop positions to check if the index is unique or an overlap
		for j in range(int(a), int(b)+1):
			if int(length_array[j]) is 0:
				length_array[j] = 1				#Use the coordinate range to change all unique indices from 0 to 1
				counter = counter + 1
			else:
				Overlapping_Index[j] = str(a)+".."+str(b)	#Dictionary (key = index, value = contig start..stop)	
		Unique_Count[str(a)+".."+str(b)] = counter			#Dictionary (key = contig start..stop, value = count of changes from 0 to 1)	
		counter = 0
	
	return Overlapping_Index, Unique_Count, length_array	

####################################################################################
#	Checks for gaps in the reference.
#	Utilizes the array of ones created from fasta file. [ 1=GAP 0=contig coverage]
####################################################################################				
def findGaps(length_array):			
	Gap_Index = {}
	for item in range(len(length_array)):
		if int(length_array[item]) is 0:
			Gap_Index[item] = "gap"				#Dictionary (key = index, value = "gap")
	return Gap_Index, length_array

########################################################################################
#	Combines the indices of overlap and gap regions into one dictionary (Unused_Index)	
########################################################################################	
def combineOverlapAndGap(Overlapping_Index, Gap_Index):			
	Unused_Index = {}			
	Unused_Index = Overlapping_Index.copy()
	Unused_Index.update(Gap_Index)
	return Unused_Index

###########################################################################################
#	Basic counts and summations of the total number of unique contigs and gap information 
###########################################################################################	
def contigCount(Unique_Count, Gap_Index):
	sum = int(0)					#Adds contigs (longest to shortest) to find total length 
	count = int(0)					#Counts total # of contigs (does not include gaps)
	total_contig_length = int(0)	#Counts total # of contigs and gaps
	for k,v in sorted(Unique_Count.items(), key=lambda k_v: (k_v[1],k_v[0]), reverse=True):
		sum = int(sum) + int(v)
		count = count + 1
	gap_count = len(Gap_Index)
	total_contig_length = sum
	return sum, count, gap_count, total_contig_length		

####################################################################################
#	Finds the median of an integer
####################################################################################	
def findMedian(length):
	return int(round(length * 0.5))

######################################################################################
#	Calculates U50/UG50 from a median value and a dictionary input 
#	Dictionary example: {startPos..stopPos : # of unique bp covered by contig}, 
#					    {0..1000 : 1000} or {500..1000 : 200}
######################################################################################	
def stats(median, Unique_Count):
	sum = int(0)
	count = int(0)					
	last = len(Unique_Count)
	for k,v in sorted(Unique_Count.items(), key=lambda k_v: (k_v[1],k_v[0]), reverse=True):
		sum = int(sum) + int(v)
		count = count + 1
		if sum >= median:
			contig = k
			value = v
			break
		if count == last and sum < median:			#Running sum of contigs never reaches median
			contig = 0
			value = 0
			sum = 0
			count = 0
	return contig, value, sum, count
	
#--------------------------------------------[ Main ]------------------------------------------------------#

fasta_length, length_array = getfasta(fasta_file)						
start_coordinates, stop_coordinates, length_list, contig_dict = getBED(BED_file)

#All lists are sorted by the length_list values in ascending order
sorted_length = [x for (x,y,z) in sorted(zip(length_list, start_coordinates, stop_coordinates))]
sorted_start = [y for (x,y,z) in sorted(zip(length_list, start_coordinates, stop_coordinates))]
sorted_stop = [z for (x,y,z) in sorted(zip(length_list, start_coordinates, stop_coordinates))]

#All lists are sorted by the length_list values in descending order
reverseList(sorted_length)
reverseList(sorted_start)
reverseList(sorted_stop)

#Find index of overlapping reads
Overlapping_Index, Unique_Count, length_array = findOverlaps(sorted_start, sorted_stop, length_array)

#Find index of gaps
Gap_Index, length_array = findGaps(length_array)

#Combine Overlapping_Index and Gap_Index into one dictionary
Unused_Index = combineOverlapAndGap(Overlapping_Index, Gap_Index)

#Return counts
contig_length, contig_count, gap_count, total_contig_length = contigCount(Unique_Count, Gap_Index)

#Calculate the median value
ref_median = findMedian(fasta_length)
contig_U50_median = findMedian(total_contig_length)

#Calculate U50/L50
contig1_from_contig, value1_from_contig, sum1_from_contig, count1_from_contig = stats(contig_U50_median, Unique_Count)	#Pass in median integer and dictionary
#Calculate UG50/LG50
contig1_from_ref, value1_from_ref, sum1_from_ref, count1_from_ref = stats(ref_median, Unique_Count)


#-------------------------------------------[ Printing Stats ]-----------------------------------------------------#

print(  "\n\nATTENTION: All coordinates are based on a ZERO-based counting system!" + \
		"\n\nReference file: " + fasta_file + \
		"\nSorted BED file: " + BED_file + \
		
		"\n\n\n\t\t\tREFERENCE" + \
		"\n\nTotal reference length: " + str(fasta_length) + \
		"\nReference median: " + str(ref_median) + \
		"\nTotal number of gaps: " + str(gap_count) + \
		
		"\n\n\n\t\t\tUNIQUE" + \
		"\n\nTotal number of unique contigs: " + str(contig_count) + \
		"\nTotal bases covered by unique contigs: " + str(contig_length) + \
		"\nUnique contig median: " + str(contig_U50_median) + \

		"\n\nU50  (Contig length): " + str(value1_from_contig) + \
		"\n      Contig coordinates: " + str(contig1_from_contig) + \
		"\n      Total length of all contigs >= U50: " + str(sum1_from_contig) + \
		"\nUL50 (Total contig count): " + str(count1_from_contig) + \
		
		"\n\nUG50  (Contig length): " + str(value1_from_ref) + \
		"\n       Contig coordinates: " + str(contig1_from_ref) + \
		"\n       Total length of all contigs >= UG50: " + str(sum1_from_ref) + \
		"\nULG50 (Total contig count): " + str(count1_from_ref) + \
		"\nUG50%: " + str('{percent:.2%}'.format(percent=(float(value1_from_ref)/float(fasta_length)))) + \
		
		"\n\n\n\t\t\tSUMMARY" + \
		"\n\nU50: \t" + str(value1_from_contig) + "\t\tUG50: \t" + str(value1_from_ref) + \
		"\nUL50: \t" + str(count1_from_contig) + "\t\tULG50: \t" + str(count1_from_ref) + \
		"\nUG50%: " + str('{percent:.2%}'.format(percent=(float(value1_from_ref)/float(fasta_length))))
		)
		
		
#-------------------------------------------[ Writing Output Files ]-----------------------------------------------------#

with open(BED_file.split(".")[0]+"_Assembly_Statistics.txt", "w") as outfile:
	outfile.write(  "\n\nATTENTION: All coordinates are based on a ZERO-based counting system!" + \
					"\n\nReference file: " + fasta_file + \
					"\nSorted BED file: " + BED_file + \
					"\n\n\n\t\t\tREFERENCE" + \
					"\n\nTotal reference length: " + str(fasta_length) + \
					"\nReference median: " + str(ref_median) + \
					"\nTotal number of gaps: " + str(gap_count) + \
					"\n\n\n\t\t\tUNIQUE" + \
					"\n\nTotal number of unique contigs: " + str(contig_count) + \
					"\nTotal bases covered by unique contigs: " + str(contig_length) + \
					"\nUnique contig median: " + str(contig_U50_median) + \
					"\n\nU50  (Contig length): " + str(value1_from_contig) + \
					"\n      Contig coordinates: " + str(contig1_from_contig) + \
					"\n      Total length of all contigs >= U50: " + str(sum1_from_contig) + \
					"\nUL50 (Total contig count): " + str(count1_from_contig) + \
					"\n\nUG50  (Contig length): " + str(value1_from_ref) + \
					"\n       Contig coordinates: " + str(contig1_from_ref) + \
					"\n       Total length of all contigs >= UG50: " + str(sum1_from_ref) + \
					"\nULG50 (Total contig count): " + str(count1_from_ref) + \
					"\nUG50%: " + str('{percent:.2%}'.format(percent=(float(value1_from_ref)/float(fasta_length))))
					)

with open(BED_file.split(".")[0]+"_gaps.txt", "w") as gap_outfile:
	gap_outfile.write("#**** ATTENTION: All coordinates are based on the BED file's ZERO-based counting system ****")
	gap_outfile.write("\n\nIndex of gap : Type\n")
	gap_outfile.writelines('{} : {}\n'.format(k,v) for k,v in sorted(Gap_Index.items()))		

with open(BED_file.split(".")[0]+"_overlaps.txt", "w") as overlap_outfile:
	overlap_outfile.write("#**** ATTENTION: All coordinates are based on the BED file's ZERO-based counting system ****")
	overlap_outfile.write("\n\nIndex of overlap : Contig coordinates\n")
	overlap_outfile.writelines('{} : {}\n'.format(k,v) for k,v in sorted(Overlapping_Index.items()))

with open(BED_file.split(".")[0]+"_contigs.txt", "w") as contig_outfile:
	contig_outfile.write("#**** ATTENTION: All coordinates are based on the BED file's ZERO-based counting system ****")
	contig_outfile.write("\n\nContig coordinates : length (only unique bp coverage, no overlaps were counted)\n")
	contig_outfile.writelines('{} : {}\n'.format(k,v) for k,v in sorted(Unique_Count.items(), key=lambda k_v: (k_v[1],k_v[0]), reverse=True))
