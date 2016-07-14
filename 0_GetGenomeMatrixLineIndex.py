# IMPORT
import sys

# MAIN

print'''

Reads a GenomeMatrix file [1] and creates a 
index of character offsets for each line
in the file and the position of the line

Optionaly, the root name of outfile can be specified 
[2, default is the infile]

IMPORTANT: This will probably take a LONG
time to run, but will speed up previous
reading of the file by using the index to
seek to specific lines

IMPORTANT: Position will stored sequentially
in the index such that position = line index + 1
The position itself is not sotred in the index
to save space

See: http://stackoverflow.com/questions/620367/how-to-jump-to-a-particular-line-in-a-huge-text-file

'''

# Read Arguments
file_lines = [ln for ln in open(sys.argv[1],"r").readlines()]

outfile = sys.argv[1]
if len(sys.argv) > 2:
    outfile = sys.argv[2]

# Get line offsets for each position
line_offsets = {}
key_list = []
offset = 0
for ln in file_lines:
    [contig,position] = ln.split("\t")[0:2]
    if contig in line_offsets.keys():
         line_offsets[contig].append(offset) 
    else:
         line_offsets[contig] = [offset]
    offset += len(ln)

# For each contig, write a different file
contig_files = {}
for contig in line_offsets:
    #Get outline 
    outlines = []
    for v in line_offsets[contig]:
        outlines.append(str(v) + "\n")

    # Write outfile
    o_file = outfile + ".index." + contig
    output = open(o_file,"w")    
    output.write("".join(outlines))
    output.close()

    # Record outfile
    contig_files[contig] = o_file


# Write a ".ref" coordinating contig and contig index
ref_outlines = []
ref_outlines.append("Target\t" + sys.argv[1] + "\n")
for contig in contig_files.keys():
    ref_outlines.append(contig + "\t" + contig_files[contig] + "\n")
output = open(outfile + ".index.ref","w")
output.write("".join(ref_outlines))
output.close()
