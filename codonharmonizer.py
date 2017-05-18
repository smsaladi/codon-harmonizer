#!/usr/bin/python3.5

import sys,getopt
import argparse
import os
import string, random
import zipfile

gencode_1 = {
"TTT":"F","TCT":"S","TAT":"Y","TGT":"C",
"TTC":"F","TCC":"S","TAC":"Y","TGC":"C",
"TTA":"L","TCA":"S","TAA":"-","TGA":"-",
"TTG":"L","TCG":"S","TAG":"-","TGG":"W",
"CTT":"L","CCT":"P","CAT":"H","CGT":"R",
"CTC":"L","CCC":"P","CAC":"H","CGC":"R",
"CTA":"L","CCA":"P","CAA":"Q","CGA":"R",
"CTG":"L","CCG":"P","CAG":"Q","CGG":"R",
"ATT":"I","ACT":"T","AAT":"N","AGT":"S",
"ATC":"I","ACC":"T","AAC":"N","AGC":"S",
"ATA":"I","ACA":"T","AAA":"K","AGA":"R",
"ATG":"M","ACG":"T","AAG":"K","AGG":"R",
"GTT":"V","GCT":"A","GAT":"D","GGT":"G",
"GTC":"V","GCC":"A","GAC":"D","GGC":"G",
"GTA":"V","GCA":"A","GAA":"E","GGA":"G",
"GTG":"V","GCG":"A","GAG":"E","GGG":"G" }

AAs = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"]

def parse_options():
	usage = "\ncodonharm.py -f <(multi)fasta_file> -o <output file> -s <frequency_file> -t <frequency_file>,<frequency_file> etc.."
	parser = argparse.ArgumentParser(usage=usage, description='Harmonize your genes for a target organism. See codonfrequencies_from_cds.py to generate frequency files.')

	input_group = parser.add_argument_group('Input')
	input_group.add_argument("-f","--fasta",  dest="fasta_filepath", help="DNA (multi)fasta file", required=True, metavar="FASTA")
	input_group.add_argument("-s","--source", dest="source",         help="Source Organism eg. Eco_MG1655")
	input_group.add_argument("-t","--target", dest="targets",        help="Target organism(s) eg. Eco_MG1655. Can be a comma separated list.", required=True, action="append")
	input_group.add_argument("-o","--output", dest="output_file",    help="Output filename (.zip)", required=True,  metavar="NAME")

	inputs = parser.parse_args()


	return inputs

def get_fasta_sequences(fasta_contents):
	"""Given a DNA fasta file, return a dictionary of the entries"""
	valid_bases = "ATCG"
	try:
		sequence_list = {}
		header = ""
		for line in fasta_contents:
			if line[0] == '>':
				header = line.strip()[1:].split(" ")[0]
				sequence_list[header] = ""
			else:
				clean_seq = line.strip().upper().replace("U","T")

				# Check for valid DNA Characters and if header is present
				if not all(char in valid_bases for char in clean_seq) or header=="":
					raise ValueError()

				sequence_list[header] += clean_seq

		#Check if there are actually sequences.
		for header in sequence_list:
			if not sequence_list[header]:
				raise ValueError()

		return sequence_list

	except ValueError as err:
		print("Not a valid DNA fasta file, maybe contains non DNA letters")
		sys.exit()

def translate(sequence,header):
	""" Return the translated protein from a sequence reading frame +1 """
	codon_sequence = split_to_codons(sequence,header)
	translated_seq = ""
	for codon in codon_sequence:
		translated_seq += gencode_1[codon]
	return translated_seq

def get_codon_freq(freq_file_path):
	""" returns relative frequencies and counts of all the codons from the given organism in a dictionairy"""
	try:
		freq_file = open(freq_file_path).readlines()
	except (OSError, IOError) as err:
		print("Invalid frequency file: ",freq_file_path)
		sys.exit()

	## Dictionairy build: counts["G"] = [['GGA', 'GGC', 'GGG'], [10216, 39284, 14444]]
	counts = {}
	name="no_name"
	for linenr, line in enumerate(freq_file):
		#FIRST GET NAME IF THERE (on 1st line: NAME,example_name)
		if linenr==0 and len(line.strip().split(","))==2:
			name = line.strip().split(",")[1]
		else:
			freq_line = line.strip().split(",")

			codon = freq_line[0]
			AA = freq_line[1]
			CodonCount = int(freq_line[2])

			if AA in counts:
				counts[AA][0].append(codon)
				counts[AA][1].append(CodonCount)
			else:
				counts[AA] = [[codon],[CodonCount]]

	frequencies = {}
	for AA in counts:
		MaxCodonCount = max(counts[AA][1])
		frequencies[AA] = [counts[AA][0],[]]
		for CodonCount in counts[AA][1]:
			CodonFreq = CodonCount/MaxCodonCount
			frequencies[AA][1].append(CodonFreq)
	return [name, frequencies, counts]

def split_to_codons(sequence,header):
	"""returns a codon list from a sequence reading frame +1 """
	if len(sequence) % 3 == 0:
		codon_sequence_list = [sequence[i:i+3] for i in range(0,len(sequence),3)]
		return codon_sequence_list
	else:
		print("Partial sequence >"+header+" not divisible by complete codons")
		sys.exit()


def get_input_sequence_freq(fasta_contents,freqs):
	""" returns relative frequencies of the input gene sequence only """
	input_sequences = get_fasta_sequences(fasta_contents)
	seq_freqs = {}
	for header in input_sequences:
		seq_freqs[header] = []

		codon_sequence = split_to_codons(input_sequences[header],header)
		AA_seq = translate(input_sequences[header],header)

		for seq_pos in range(len(AA_seq)):
			codon_pos = [i for i,x in enumerate(freqs[AA_seq[seq_pos]][0]) if x == codon_sequence[seq_pos]][0]
			seq_freqs[header].append([codon_sequence[seq_pos],AA_seq[seq_pos],freqs[AA_seq[seq_pos]][1][codon_pos]])

	return seq_freqs


def harmonize_sequences(source_gene_freqs,target_freqs):
	#TODO: Add Restriction sites.
	harmonized_sequence_freqs = []

	for AA_position in source_gene_freqs:
		AA = AA_position[1]
		source_codon_freq = AA_position[2]

		target_AA_codon  = target_freqs[AA][0]
		target_AA_freqs = target_freqs[AA][1]

		closest_freq = min(enumerate(target_AA_freqs), key=lambda x: abs(x[1]-source_codon_freq))[0]

		harmonized_sequence_freqs.append([target_freqs[AA][0][closest_freq],AA,target_freqs[AA][1][closest_freq]])

	return harmonized_sequence_freqs

def get_index_score(freqs,native_freq):
	''' calculate index score or CHI, Codon Harmonization Index '''
	seq_len = len(native_freq)
	diffsum=0

	for AApos in range(seq_len):
		diffsum+=abs(freqs[AApos][2]-native_freq[AApos][2])

	score = 1/seq_len*diffsum
	return score

def CAI(gene_freqs):
	'''calculate average frequencies '''
	gene_length = len(gene_freqs)
	total = 0
	for f in gene_freqs: total+=f[2]
	cai = total/gene_length
	return cai


def main(argv):
	inputs = parse_options()

	# READ FASTA FILE
	try: fastafile = open(inputs.fasta_filepath,"r").readlines()
	except (OSError, IOError) as err:
		print("Unable to open input fasta file")
		sys.exit()

	# GET SOURCE FREQUENCIES
	source_codon_freqs = get_codon_freq(inputs.source)
	source_name, source_freqs, source_counts = source_codon_freqs[0], source_codon_freqs[1], source_codon_freqs[2]
	source_gene_freqs = get_input_sequence_freq(fastafile, source_freqs)

	# GET TARGET FREQUENCIES
	target_gene_freqs = {}
	target_freqs = {}
	target_counts = {}

	input_targets = []
	input_targets += sum([i.split(",") for i in inputs.targets],[])

	for target in input_targets:
		target_codon_freqs = get_codon_freq(target)
		target_name,target_freq,target_count = target_codon_freqs[0], target_codon_freqs[1], target_codon_freqs[2]

		target_freqs[target_name] = target_freq
		target_counts[target_name] = target_count
		target_gene_freqs[target_name] = get_input_sequence_freq(fastafile,target_freqs[target_name])


	output_zip = zipfile.ZipFile(inputs.output_file, "w")

	# LOOP THROUGH INPUT SEQUENCES
	for sequence in source_gene_freqs:
		for target in target_gene_freqs:
			harm_out_filename = sequence+"_"+source_name+"-"+target+"_harmonized.csv"
			harm_out_file = open(harm_out_filename,"w")
			harm_gene_freqs = harmonize_sequences(source_gene_freqs[sequence],target_freqs[target])

			source_sequence     = "".join([AA_pos[0] for AA_pos in source_gene_freqs[sequence]])
			harmonized_sequence = "".join([AA_pos[0] for AA_pos in harm_gene_freqs])

			CAI_native = CAI(source_gene_freqs[sequence])
			CAI_nonNative = CAI(target_gene_freqs[target][sequence])
			CAI_harmonized = CAI(harm_gene_freqs)

			CHI = get_index_score(harm_gene_freqs,source_gene_freqs[sequence])
			non_CHI = get_index_score(target_gene_freqs[target][sequence],source_gene_freqs[sequence])

			# WRITE TO CSV FILE
			harm_out_file.write("Source Organism:,"+source_name+"\n")
			harm_out_file.write("Target Organism:,"+target+"\n")
			harm_out_file.write("Gene name:,"+sequence+"\n")
			harm_out_file.write("Source Gene:,"+source_sequence+"\n\n")

			harm_out_file.write("Harmonized gene:,"+harmonized_sequence+"\n\n")

			harm_out_file.write(",,,,,CHI:,"+str(non_CHI)+",,CHI:,"+str(CHI)+"\n")
			harm_out_file.write(",,CAI:,"+str(CAI_native)+",,CAI:,"+str(CAI_nonNative)+",,CAI:,"+str(CAI_harmonized)+"\n\n")

			harm_out_file.write(","+source_name+" Source gene in "+source_name+",,,"+source_name+" Source gene in "+target+",,,"+source_name+" harmonized gene in "+target+"\n")
			harm_out_file.write("AA pos.,codon,AA,codon freq.,codon,AA,codon freq.,codon,AA,codon freq.\n")

			for i,source_AA_pos in enumerate(source_gene_freqs[sequence]):
				# native frequencies
				harm_out_file.write(str(i+1)+",")
				harm_out_file.write(source_AA_pos[0]+",")
				harm_out_file.write(source_AA_pos[1]+",")
				harm_out_file.write(str(source_AA_pos[2])+",")

				# non native frequencies
				harm_out_file.write(target_gene_freqs[target][sequence][i][0]+",")
				harm_out_file.write(target_gene_freqs[target][sequence][i][1]+",")
				harm_out_file.write(str(target_gene_freqs[target][sequence][i][2])+",")

				# harmonized frequencies
				harm_out_file.write(harm_gene_freqs[i][0]+",")
				harm_out_file.write(harm_gene_freqs[i][1]+",")
				harm_out_file.write(str(harm_gene_freqs[i][2])+"\n")


			# WRITE ALL CODON FREQUENCIES OF SOURCE AND TARGETS
			harm_out_file.write("\n\n,,"+source_name+",,"+target+"\n")
			harm_out_file.write("AA,codon,source organism freq.,source organism total,target organism freq,target organism total.\n")

			for AA in AAs:
				for i,codon in enumerate(source_counts[AA][0]):
					harm_out_file.write(AA+","+codon+","+str(source_freqs[AA][1][i])+","+str(source_counts[AA][1][i])+","+str(target_freqs[target][AA][1][i])+","+str(target_counts[target][AA][1][i])+"\n")


			harm_out_file.close()
			output_zip.write(harm_out_filename)
	output_zip.close()

if __name__ == "__main__":
	main(sys.argv[1:])
