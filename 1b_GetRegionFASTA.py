# IMPORT
import sys

# MAIN

print'''

Takes a list regions [1] and a *.index.ref [2] and returns
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

# Convert the regional information into sequences and write
# Note: This requires tranposing rows and columns
for region in reg_dict:

    # Get lines and convert to sequence
    region_lines = [ln.strip().split("\t")[2:] for ln in reg_dict[region]]
    fasta_lines = []
    for index in range(len(region_lines[0])):
        fasta_lines.append("> " + str(index)+"\n")
        fasta_lines.append("".join(ln[index] for ln in region_lines)+"\n")

    # Write fasta file
    output = open(sys.argv[1] + "." + "_".join(list(region)) + ".fasta","w")
    output.write("".join(fasta_lines))
    output.close()
    
