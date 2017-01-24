# GenomeMatrixProcessing
Python functions and scripts for parsing regions out of a genome matrix file (see  http://1001genomes.org/data/MPI/MPICao2010/releases/current/genome_matrix/) and converting
to other another format (i.e. FASTA) or doing analysis directly on the extracted region

# Useage 
The central script in this pipeline is "0_GetGenomeMatrixLineIndex.py" which for each sequence in the genome matrix, created a single column index [target_file].index.[seq] 
where there value of the Nth line is the offset for reading the line in the genome matrix file which represent the Nth base of a chromsome. A reference file, 
[target_file].index.ref while recording the original target file and the relationship between sequence and each index

Once the index for the genome matrix file has been created, the following script can be used to process the genome matrix:
  1. 1a_GetRegions.py: extracts regions from the the genome matrix file into a single file where each region is seperated by a header line beginning with ">"
     *The regions file is a three column file with the following format: [seq]\t[start]\t[stop]\n
     *The name of each region in the header line will be [seq]_[start]_[stop]
  2. 1b_GetRegionFASTA.py: extracts the base sequeces of region from the genome matrix and write the sequences of each region into a seperate multi-FASTA file
     *Uses the same regions file format as above
     *Outfile files have the following format [region_file].[seq]_[start]_[stop].fasta
     *Sequences within each fasta are named sequencetial, beginning with "0" as the referrence sequence
  3. 1c_GetRegionStats.py: calculates stats (Nt Diversity and Tajima's D) for the sequences in a defined region
     *Uses the same regions file format as above
     *Outfiles have the following format: [Region]\t[NtDiversity]\t[TajimasD]\n
     *If there are no differential sites in a region, "NoDiffSites" will in that cell instead of a float value
  4. 1d_GetRegionStats_fromBED.py: preforms the same function as 1c_GetRegionStats.py but uses a standard BED format file instead of a regions file
      *BED format file should be in the following format: [Seq]\t[Start]\t[Stop]\t[SeqName]\n. All four values will be used to name regions
      *Output format is the same as 1c_GetRegionStats.py
