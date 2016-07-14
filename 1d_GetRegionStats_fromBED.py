# IMPORT
import sys

# FUNCTIONS

def SeqPairDiff(seq1,seq2):
    ''' Finds the aligned length, location of differences, and number of differences
    between two genomic sequences

    '''
    bases = ("A","C","G","T")
    length = 0
    sites = []
    difference = 0

    for index in range(len(seq1)):
        [base1,base2] = [seq1[index],seq2[index]]
        if base1 in bases and base2 in bases:
            length += 1
            if not base1 == base2:
                difference += 1   
                sites.append(index)

    return length, sites, difference

def RegionToSeq(region_lines):
    ''' Convert a list of bases for a region from a GenomeMatrix file into
    a lsit of sequences

    '''

    sequences = []
    trunc_lines = [ln.strip().split("\t")[2:] for ln in region_lines]
    for index in range(len(trunc_lines[0])):
	sequences.append([ln[index] for ln in trunc_lines])

    return sequences

def DifferenceVaraince(n,S):
    ''' Calculates the variance in the expected number of SNPs between
    sequences for Tajima's D where "n" in the number of sequences and
    S in the number of segregating sites. 

    See https://en.wikipedia.org/wiki/Tajima%27s_D

    '''

    import math

    a1 = sum([1/float(i+1) for i in range(n-1)])
    a2 = sum([1/float((i+1)**2) for i in range(n-1)])

    b1 = float(n+1)/float(3*(n-1))
    b2 = float(2*(n**2 + n + 3))/float(9*n*(n-1))
    
    c1 = b1 - (1/float(a1))
    c2 = b2 - float(n+2)/float(a1*n) + float(a2)/float(a1**2)
 
    e1 = c1/float(a1)
    e2 = c2/float(a1**2 + a2)

    V_of_d = math.sqrt(e1*S + e2*S*(S-1))

    return V_of_d 

def SeqStats(seqs):
    ''' Takes a list of sequences and calculates stats for these sequences

        Current does NucleotideDiversity and TajimasD
    '''

    # Get the number of sequences
    seq_num = len(seqs)

    # Make list to store our stats
    seq_lens = []  # List of the aligned length between sequences
    seq_sites = [] # List of sites with different SNPs
    seq_diff =  [] # List of differences between sequences
    
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            [lens,sites,diffs] = SeqPairDiff(seqs[i],seqs[j])
            seq_lens.append(lens)
            seq_sites.extend(sites)
            seq_diff.append(diffs)

    seq_sites = list(set(seq_sites))

    #Test code
    #print "# of Diff Sites: " + str(len(seq_sites))
    #print "Seq Diffs: " + str(seq_diff)

    # Calculate Stats #

    if len(seq_sites) > 0:
        # Nucleotide Diversity: See https://en.wikipedia.org/wiki/Nucleotide_diversity
        seq_frequency = 1/float(seq_num)
        sum_diffs = sum([float(seq_diff[i])/float(seq_lens[i]) for i in range(len(seq_lens)) if seq_lens[i] > 0])
        nt_diversity = 2 * (seq_frequency**2) * float(sum_diffs)     

        #Test Code
        #print "SeqFrequency: " + str(seq_frequency)
        #print "Sum Diffs: " + str(sum_diffs)

        # Tajima's D: https://en.wikipedia.org/wiki/Tajima%27s_D
        avg_polymorphism = float(sum(seq_diff)) / float(len(seq_diff))
        expected_polymorphism = float(len(seq_sites)) / sum([1/float(i+1) for i in range(seq_num-1)])
        stdv_of_expected = DifferenceVaraince(seq_num,len(seq_sites))
        tajimas_D = (avg_polymorphism - expected_polymorphism)/stdv_of_expected

        #Test Code
        #print "Avg Poly: " + str(avg_polymorphism)
        #print "Expected: " + str(expected_polymorphism)
    else:
        nt_diversity = "NoDiffSites"
        tajimas_D = "NoDiffSites"

    # Save values
    return [nt_diversity,tajimas_D]



# MAIN

print'''

Takes a BED file [1] and a *.index.ref [2] and returns
the relevant regions as a multi-FASTA file

'''

region_lines = [ln.strip() for ln in open(sys.argv[1],"r").readlines()]
ref_lines = [ln.strip() for ln in open(sys.argv[2],"r").readlines()]

# Make dictionary of relationship in the reference file
ref_dict = {}
for ln in ref_lines:
    [name,file] = ln.split("\t")
    ref_dict[name] = file

target_file = ref_dict["Target"]

# Make a dictioanry of regions
outlines = []
for ln in region_lines:
    [contig,start,stop,label] = ln.split("\t")
    
    # Get offset from index file
    try:
        contig_index_file = ref_dict[contig] 
    except KeyError as e:
        print "No index file for contig " + contig

 
    start_offset = open(contig_index_file,"r").readlines()[int(start)-1].strip()
    number_of_lines = int(stop) - int(start)

    # Get lines from target file
    target_read = open(target_file)    
    target_read.seek(int(start_offset))
    lines = [target_read.readline() for n in range(number_of_lines+1)]
    target_read.close()

    # Convert regions to Seq
    seqs = RegionToSeq(lines)

    # Get stats for the sequences 
    seq_stats = SeqStats(seqs)

    # Write outlines
    outlines.append("_".join([contig,start,stop,label]) + "\t" + "\t".join([str(s) for s in seq_stats]) + "\n")

# Write output 
output = open(sys.argv[1] + ".SeqStats","w")
output.write("Region\tNtDiversity\tTajimasD\n")
output.write("".join(outlines))
output.close()
