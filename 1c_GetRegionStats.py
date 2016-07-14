# IMPORT
import sys

# FUNCTIONS

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

# MAIN

print'''

Takes a list regions [1] and a *.index.ref [2] and returns
statistics for each region (NtDiversity and Tajima's D)

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
reg_dict = {}
for ln in region_lines:
    [contig,start,stop] = ln.split("\t")
    
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

    reg_dict[(contig,start,stop)] = lines

# Convert the regional information into sequences
seq_keys = []
seq_dict = {}
for region in reg_dict:

    # Get lines and convert to sequence
    sequences = []
    region_lines = [ln.strip().split("\t")[2:] for ln in reg_dict[region]]
    for index in range(len(region_lines[0])):
	sequences.append([ln[index] for ln in region_lines])

    seq_keys.append(region)
    seq_dict[region] = sequences
    
# For each region, calculate the pairwise difference between all sequenecs 
stats_dict = {}
for region in seq_keys:

    # Get sequences 
    r_seq = seq_dict[region]
    seq_num = len(r_seq)

    # Make list to store our stats
    seq_lens = []  # List of the aligned length between sequences
    seq_sites = [] # List of sites with different SNPs
    seq_diff =  [] # List of differences between sequences
    
    for i in range(len(r_seq)):
        for j in range(i+1,len(r_seq)):
            [lens,sites,diffs] = SeqPairDiff(r_seq[i],r_seq[j])
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
        variance =  DifferenceVaraince(seq_num,len(seq_sites))
        tajimas_D = (avg_polymorphism - expected_polymorphism)/variance

        #Test Code
        #print "Avg Poly: " + str(avg_polymorphism)
        #print "Expected: " + str(expected_polymorphism)
    else:
        nt_diversity = "NoDiffSites"
        tajimas_D = "NoDiffSites"

    # Save values
    stats_dict[region] = [nt_diversity,tajimas_D]

# Write output 
outlines = []
for region in seq_keys:
    [nt,D] = stats_dict[region]
    outlines.append("_".join(list(region)) + "\t" + str(nt) + "\t" + str(D) + "\n")
output = open(sys.argv[1] + ".SeqStats","w")
output.write("Region\tNtDiversity\tTajimasD\n")
output.write("".join(outlines))
output.close()
