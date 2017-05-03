import sys
import os

# refersence files:
germline_birds_file = '/home/users/a.steep/databases/samples/germline_sample_rnaseq_list_017NNN-N_SN.txt'
tumor_birds_file = '/home/users/a.steep/databases/samples/tumor_sample_rnaseq_list_017NNN-N.txt'

# Gene ID to gene name file
gene_ann_file = '/home/proj/MDW_genomics/steepale/gene_fusions/data/gal5.gene_id.gene_name.tsv'

# build a dictionary for germline output files
# Note: only germline samples of final filtered calls provided
germline_results = {}
for gbird in open(germline_birds_file):
	gbird = gbird.rstrip()
	germline_results[gbird] = '/home/proj/MDW_genomics/steepale/gene_fusions/data/chimericJunctions_' + gbird + '.txt'

# build a dictionary for tumor output files
tumor_results = {}
for tbird in open(tumor_birds_file):
	tbird = tbird.rstrip()
	tumor_results[tbird] = '/home/proj/MDW_genomics/steepale/gene_fusions/data/chimericJunction_candidates_' + tbird + '.txt'

# Build dictionary to associate geneID with gene name
gene_id2gene_name = {}
for gene_line in open(gene_ann_file):
	gene_line = gene_line.rstrip()
	geneID = gene_line.split('\t')[0]
	gene_name = gene_line.split('\t')[1]
	gene_id2gene_name[geneID] = gene_name

# Outfile
outfile = open('/home/proj/MDW_genomics/steepale/gene_fusions/data/chimeric_junctions_unfiltered_annotated.int', 'w')

# Write header to outfile
outfile.write('#juncCoord' + '\t' + 'gene_name_A' + '\t' + 'gene_name_B' + '\t' + 'sample_num' + '\t' + 'samples' + '\t' + 'totalNbPE' + '\t' + 'nbSpanningPE' + '\t' + 'nbDiscordantPE' + '\n')

# Create a dictionary for each variant and associate the number of samples with it
pos2sample = {}
# Create a dictionay of total reads across samples
pos2totalNbPE = {}
# Create a dictionary of split read number across samples
pos2split = {}
# Create a dictionary for paired read number across samples
pos2paired = {}
# Iterate over lines of tumor infiles and collect fields
for tbird, tumor_infile in tumor_results.items():
	for tline in open(tumor_infile):
		if tline.split('\t')[0] != 'juncCoord':
			tline = tline.rstrip()
			incol = tline.split('\t')
			pos = incol[0]
			totalNbPE = incol[1]
			split_read_num = incol[7]
			paired_read_num = incol[2]
			sample = str(tbird.rstrip())
			# Append the list in the junction to sample dictionary
			if pos in pos2sample.keys():
				pos2sample[pos].append(sample)
			else:
				pos2sample[pos] = [sample]
			# Append total reads dictionary
			if pos in pos2totalNbPE.keys():
				pos2totalNbPE[pos].append(totalNbPE)
			else:
				pos2totalNbPE[pos] = [totalNbPE]
			# Append split read number dictionary
			if pos in pos2split.keys():
				pos2split[pos].append(split_read_num)
			else:
				pos2split[pos] = [split_read_num]
			# Append paired read number dictionary
			if pos in pos2paired.keys():
				pos2paired[pos].append(paired_read_num)
			else:
				pos2paired[pos] = [paired_read_num]

# Iterate over the lines of the germline files and collect fields in oder to add germline samples to samples field for later filtering
for gbird, germline_infile in germline_results.items():
	for gline in open(germline_infile):
		if gline.split('\t')[0] != 'juncCoord':
			gline = gline.rstrip()
			incol = gline.split('\t')
			pos = incol[0]
			sample = str(gbird.rstrip())
			# Append the list in the junction to sample dictionary
			if pos in pos2sample.keys():
				pos2sample[pos].append(sample)
			else:
				pos2sample[pos] = [sample]

# Now that dictionary has been constructed, iterate over tumor files again and print file

# Iterate over lines of tumor infiles and collect fields
for tbird, tumor_infile in tumor_results.items():
	for tline in open(tumor_infile):
		if tline.split('\t')[0] != 'juncCoord':
			tline = tline.rstrip()
			incol = tline.split('\t')
			pos = incol[0]
			totalNbPE = incol[1]
			split_read_num = incol[7]
			paired_read_num = incol[2]
			geneID_A = incol[23]
			if geneID_A in gene_id2gene_name.keys():
				gene_name_A = gene_id2gene_name[geneID_A]
			geneID_B = incol[24]
			if geneID_B in gene_id2gene_name.keys():
				gene_name_B = gene_id2gene_name[geneID_B]
			# Filter out any variants that are found in any of the candidate germline files
			if any(x in germline_results.keys() for x in pos2sample[pos]) == False:
				# Collect sample number
				sample_num = str(len(pos2sample[pos]))
				# format the lists in samples dictionary
				if len(pos2sample[pos]) == 1:
					samples = pos2sample[pos][0]
				elif len(pos2sample[pos]) > 1:
					samples = ";".join(pos2sample[pos])
				else:
					print('\n' + '\n' + '\n' + 'WARNING' + '\n' + '\n' + '\n')
				# format list in totalNbPE dictionary
				if len(pos2totalNbPE[pos]) == 1:
					totalNbPE_fin = pos2totalNbPE[pos][0]
				elif len(pos2totalNbPE[pos]) > 1:
					totalNbPE_fin = ";".join(pos2totalNbPE[pos])
				else:
					print('\n' + '\n' + '\n' + 'WARNING' + '\n' + '\n' + '\n')

				# format list in split read dictionary
				if len(pos2split[pos]) == 1:
					split_read_fin = pos2split[pos][0]
				elif len(pos2split[pos]) > 1:
					split_read_fin = ";".join(pos2split[pos])
				else:
					print('\n' + '\n' + '\n' + 'WARNING' + '\n' + '\n' + '\n')
				# format list in paired read dictionary
				if len(pos2paired[pos]) == 1:
					paired_read_fin = pos2paired[pos][0]
				elif len(pos2paired[pos]) > 1:
					paired_read_fin = ";".join(pos2paired[pos])
				else:
					print('\n' + '\n' + '\n' + 'WARNING' + '\n' + '\n' + '\n')
				# write to output file
				outfile.write(pos + '\t' + gene_name_A + '\t' + gene_name_B + '\t' + sample_num + '\t' + samples + '\t' + totalNbPE_fin + '\t' + paired_read_fin + '\t' + split_read_fin + '\n')
outfile.close()

