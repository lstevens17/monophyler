#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from ete3 import Tree
from Bio import SeqIO
from os.path import basename

def parse_fasta(fasta_file):
	with open(fasta_file, "rU") as fasta: # rU converts all newline characters to \n
		fasta_dict = {}
		for record in SeqIO.parse(fasta, "fasta"):
			fasta_dict[record.id] = record.seq
	return fasta_dict

def find_dups(fasta_dict):
	tmp_list = []
	dup_list = []
	for seq in fasta_dict:
		species = seq.split(".")[0]
		if not species in tmp_list:
			tmp_list.append(species)
		else:
			dup_list.append(species)
	return dup_list

def parse_tree(tree_file):
	with open(tree_file, "rU") as tree_file:
		lines = ''
		for line in tree_file:
			lines += line.rstrip("\n")
		tree = Tree(lines)
		outgroup = tree.get_midpoint_outgroup()
		tree.set_outgroup(outgroup)
		return tree

def identify_monophyly(tree, dup_list, fasta_dict, threshold):
	monophyletic_list = []
	for sp in dup_list:
		seq_list = []
		for seq in fasta_dict:
			if seq.split(".")[0] == sp:
				seq_list.append(seq)
		common_ancestor = tree.get_common_ancestor(seq_list)
		if len(common_ancestor.get_leaf_names()) == 2:
			monophyletic_list.append(seq_list[0].split(".")[0])
		elif len(common_ancestor.get_leaf_names()) == 3:
			if float(common_ancestor.support) <= threshold:
				monophyletic_list.append(seq_list[0].split(".")[0])
	return monophyletic_list

def write_file(monophyletic_list, dup_list, fasta_dict, fasta):
	if sorted(monophyletic_list) == sorted(dup_list):
		with open(fasta + "_pruned", 'w') as outfile:
			output_dict = {}
			for seq in fasta_dict:
				species = seq.split(".")[0]
				if not species in output_dict:
					output_dict[species] = seq
				else:
					existing_seq = output_dict[species]
					if len(fasta_dict[seq]) > len(existing_seq):
						output_dict[species] = seq
			for species in output_dict:
				seq = output_dict[species]
				outfile.write(">" + output_dict[species] + "\n" + str(fasta_dict[seq]) + "\n")
		print "All species monophyletic\t"
	else:
		print ",".join(set(monophyletic_list).symmetric_difference(set(dup_list))) + " not monophyletic\t"

if __name__ == "__main__":

	SCRIPT = basename(__file__)

	infile, value, outfile = '', 0,''
	try:
		tree = sys.argv[1]
		fasta = sys.argv[2]
		threshold = int(sys.argv[3])
	except:
		sys.exit("USAGE: ./%s %s %s %s" % (SCRIPT, "[tree_with_bootstraps.nwk]", "[unaligned_seqs.fa]", "[min_bootstrap_val]"))

	fasta_dict = parse_fasta(fasta)
	dup_list = find_dups(fasta_dict)
	tree = parse_tree(tree)
	monophyletic_list = identify_monophyly(tree, dup_list, fasta_dict, threshold)
	write_file(monophyletic_list, dup_list, fasta_dict, fasta)
