#!/usr/bin/env python
import argparse
def get_args():
    parser = argparse.ArgumentParser(description="Use sorted and uniquely mapped sam files. code written for single end reads and assumes a known list of UMIs is available. First read written to file if duplicates. Code does not address UMI error")
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file")
    parser.add_argument("-u", "--umi", help="designates file containing the list of UMIs")
    parser.add_argument("-o", "--outfile", help="designates absolute file path to sorted sam file")
    #parser.add_argument("-h", "--help", required=False, help="sorted and uniquely mapped sam files. code written for single end reads and assumes a known list of UMIs is available. First read written to file if duplicates. Code does not address UMI error")
    return parser.parse_args()
#print("yes i print")
def flag_revcomp(strand):
    "determine if bitwise flag in sam file indicates read is reverse complement"
    if ((strand & 16) == 16):
        rev_comp = "reverse"
    else:
        rev_comp = "forward"
    return rev_comp
import re
args=get_args()
#/Users/samobermiller/bioinfo/Bi624/Deduper-samobermiller/obermiller_deduper.py -u /Users/samobermiller/bioinfo/Bi624/Deduper-samobermiller/unit_test/STL96.txt -f /Users/samobermiller/bioinfo/Bi624/Deduper-samobermiller/unit_test/test.sam -o /Users/samobermiller/bioinfo/Bi624/Deduper-samobermiller/practice_output/practice.sam
outsam=open(args.outfile, "w")
known_umi_list=[]
#create UMI_list from provided UMI file
with open(args.umi, "r") as umi_list:
    for line in umi_list:
        #print(line)
        known_umi_list.append(line.strip())
#print(len(known_umi_list))

final_dict={}
current_spec=[]
header_count=0
uniq_reads=0
wrong_umi=0
dup_count=0
countperchrom=0
chrom_dict={}
with open(args.file, "r") as sam:
    #print("current list:")
    for line in sam:
        if line.startswith("@"):
        # if line.startswith('@HD') or line.startswith('@PG') or line.startswith('@SQ'):
            #print(line)
            header_count+=1
            outsam.write(line)
            #write header to sam file
        else:
            current_line = line
            #print(current_line)
            current_line=current_line.strip('\n')
            elements=current_line.split('\t')
            #elements=line.strip('\n')
            #print(elements)
            qname=elements[0].split(':')
            #print(qname)
            umi=qname[7]
            #print(umi)
            strand=int(elements[1])
            #print(strand)
            chromosome=elements[2]
            LM=int(elements[3])
            #print(LM)
            cigar=str(elements[5])
            #print(cigar)
            start=int(LM)
            current_spec=[umi, strand, chromosome, LM, cigar, start]
            #print(cigar)
            #current_spec = [read's UMI, strand (column2), chromosome(column3), LM_position(column4), cigar(column6), start_position=LM_position]
            if umi in known_umi_list:
                #print(umi)
                if flag_revcomp(strand) == 'forward':
                    #if strand is forward
                    #print("forward")
                    if re.match(r'^(\d+)S.*', cigar):
                        #print('left soft clipped')
                        #print(cigar)
                    #if S at begining of cigar
                        start_softclip=cigar.split('S')
                        #adjustS="##"
                        start_adjustS=int(start_softclip[0])
                        #print(start_adjustS)
                        start=start-start_adjustS
                        #print(start)
                        current_spec[5]=start
                        #print(current_spec)
                #print(current_spec)
                if flag_revcomp(strand) == 'reverse':
                #if strand is reverse
                    #print("reverse")
                    if re.match(r'.*(\d+)S$', cigar):
                    #if cigar from current_spec ends with "##s":
                        #print("right soft clipped")
                        end_softclip=re.findall(r'.*(\d+)S$', cigar)
                        #adjustS="##"
                        #print(end_softclip)
                        #print(end_softclip[0])
                        end_adjustS=int(end_softclip[0])
                        #print(end_adjustS)
                        start=end_adjustS+start
                        current_spec[5]=start
                        #print(start)
                        #start position in current_spec=adjustS+LM_position
                #print(current_spec)
                    for letter in re.findall(r'(\d+)(M{1})', cigar):
                        #if "s' not in cigar from current_spec:
                        #print('M present')
                        #print(cigar)
                        match=letter
                        #print(match)
                        match_adjust=int(match[0])
                        #print(match_adjust)
                        start+=match_adjust
                        #print(start)
                        current_spec[5]=start
                        #print(current_spec)
                        #start_position in current_spec = adjustM + LM_position
                    #print(current_spec)
                    for letter in re.findall(r'(\d+)(N{1})', cigar):
                        #if cigar from current_spec has "N":
                        #print('N present')
                        #print(cigar)
                        N_value=letter
                        #print(N_value)
                        adjust_N=int(N_value[0])
                        #print(adjust_N)
                        start=adjust_N+start
                        #print(start)
                        current_spec[5]=start
                        #start_position in current_spec = adjustN + LM_position
                        #print(current_spec)
                    for letter in re.findall(r'(\d+)(D{1})', cigar):
                        #print("D present")
                        #print(cigar)
                        D_value=letter
                        adjust_D=int(D_value[0])
                        #print(adjust_D)
                        start=adjust_D+start
                        #print(start)
                        current_spec[5]=start
                        #print(current_spec)
                        #start_position in current_spec = adjustD + LM_position
                #print(current_spec[5])
                #current_spec[5]=current_spec[5]-1
                current_spec=(current_spec[0],current_spec[1],current_spec[2],current_spec[5])
                #print(current_spec)
                if current_spec not in final_dict:
                    uniq_reads+=1
                    final_dict[current_spec] = 1
                    #print("new umi in final dict")
                    #print(final_list)
                    if current_spec[2] in chrom_dict:
                        chrom_dict[current_spec[2]]+=1
                    else:
                        chrom_dict[current_spec[2]]=1
                    outsam.write(current_line+'\n')
                    current_spec=[]
                    continue
                else:
                    #print("duplicate found")
                    final_dict[current_spec] += 1
                    dup_count+=1
                    current_spec=[]
                    continue
            else:
                wrong_umi+=1
                #print("umi not in known list")
                        #print("yes")
with open("./requirements.txt", "w") as r:
    r.write("# header lines:"+'\n')
    r.write(str(header_count)+'\n')
    r.write("# uniq reads:"+'\n')
    r.write(str(uniq_reads)+'\n')
    r.write("# wrong UMIs"+'\n')
    r.write(str(wrong_umi)+'\n')
    r.write("# dups removed"+'\n')
    r.write(str(dup_count)+'\n')
    r.write("reads per chrom:"+'\n')
    for key, value in chrom_dict.items():
        count=f'{key}\t{value}\n'
        r.write(str(count))