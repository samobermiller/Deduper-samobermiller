#!/usr/bin/env python
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="a program for kmer")
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file")
    parser.add_argument("-u", "--umi", help="designates file containing the list of UMIs")
    parser.add_argument("-o", "--outfile", help="designates absolute file path to sorted sam file")
    #add -h help
    #parser.add_argument("-h", "--help", help="prints a USEFUL help message (see argparse docs)")
    return parser.parse_args()
args=get_args()
#/Users/samobermiller/bioinfo/Bi624/Deduper-samobermiller/obermiller_deduper.py -u /Users/samobermiller/bioinfo/Bi624/Deduper-samobermiller/unit_test/STL96.txt -f /Users/samobermiller/bioinfo/Bi624/Deduper-samobermiller/unit_test/test.sam -o /Users/samobermiller/bioinfo/Bi624/Deduper-samobermiller/practice_output/practice.sam
outsam=open(args.outfile, "w")
known_umi_list=[]
#create UMI_list from provided UMI file
with open(args.umi, "r") as umi_list:
    for line in umi_list:
        #print(line)
        known_umi_list.append(line.strip())
# create empty dictionary (final_list): {()} will hold UMI, position, chrom, \
# strand for each read and will be used to compare future reads for duplicates, \
# don't need value
final_list={()}
# create temporary array (current_spec): will hold current UMI, position, chrom, strand
current_spec=[]
with open(args.file, "r") as sam:
    for line in sam:
        if line.startswith('@HD') or line.startswith('@PG') or line.startswith('@SQ'):
            #print(line)
            outsam.write(line)
            #write header to sam file
        else:
            if line=="":
                break
            else:
                current_line = line
                elements=line.split('\t')
                #print(elements)
                qname=elements[0].split(':')
                #print(qname)
                umi=qname[7]
                #print(umi)
                strand=elements[1]
                #print(strand)
                chromosome=elements[2]
                LM=elements[3]
                cigar=elements[5]
                start=LM
                current_spec=[umi, strand, chromosome, LM, cigar, start]
                #print(current_spec)
#             * current_spec = [read's UMI, strand (column2), chromosome(column3), LM_position(column4), cigar(column6), start_position=LM_position]
#             * if UMI in current_spec is in UMI_list:
#                 * if strand_function applied to strand in current_spec = mapped:
#                     * if bitwise_revcomp = false:
#                         * if cigar from current_spec starts with "##s":
#                             * adjustS="##"
#                             * start_position in current_spec=LM_position-adjustS
#                         * if cigar from current_spec ends with "##s":
#                             * continue
#                         * if "s' not in cigar from current_spec:
#                             * continue
#                     * if bitwise_revcomp applied to strand in current_spec = true:
#                         * if cigar from current_spec starts with "##s":
#                             * continue
#                         * if cigar from current_spec ends with "##s":
#                             * adjustS="##"
#                             * start position in current_spec=adjustS+LM_position
#                         * if "s' not in cigar from current_spec:
#                             * adjustM= ## before M
#                             * start_position in current_spec = adjustM + LM_position
#                         * if cigar from current_spec has "N":
#                             * adjustN= ## before N
#                             * start_position in current_spec = adjustN + LM_position
#                         * if cigar from current_spec has "D":
#                             * adjustD= ## before D
#                             * start_position in current_spec = adjustD + LM_position
#                 * is UMI, start_position, chrom, strand in final array?
#                     * if yes:
#                         * empty current_line and current_spec
#                         * break
#                     * if no:
#                         * add info from current_spec to final array
#                         * write entire line from current_line to final sam file
#                         * empty current_line and current_spec
#                         * break
#                 * if strand_function == unmapped:
#                     * empty current_line and current_spec
#                     * break
#             * if UMI in current_spec array not in UMI_list:
#                 * empty current_line and current_spec
#                 * break

# - Determine high level functions

# def strand_function:
#     "determine if bitwise flag in sam file indicates read is mapped (true) or unmapped (false)"
#         return mapped
# assert strand_function("4")==unmapped

# def bitwise_revcomp:
#     "determine if bitwise flag in sam file indicates read is reverse complement"
#     return true
# assert bitwise_revcom("16")==true