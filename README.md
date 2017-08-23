# Monophyler.py 

Reimplementation of some features in PhyloTreePruner, as it's problematic to run and has some annoying input requirements. 

Given a tree with two loci for one or more species, it will determine if both loci are monophyletic/in-paralagous and spit out the largest sequence of the two, yielding a FASTA of single-copy orthologs

USAGE:
	`./monophyler.py [tree_with_bootstraps.nwk] [unaligned_seqs.fa] [min_bootstrap_val]`

EXAMPLE:
	`./monophyler.py test_data/RAxML_bipartitions.OG0004442.fa.aln test_data/OG0004442.fa 50`



