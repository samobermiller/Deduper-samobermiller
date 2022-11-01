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
import bioinfo
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

final_list={}

# create temporary array (current_spec): will hold current UMI, position, chrom, strand
current_spec=[]
with open(args.file, "r") as sam:
    for line in sam:
        if line.startswith('@HD') or line.startswith('@PG') or line.startswith('@SQ'):
            #print(line)
            outsam.write(line)
            #write header to sam file
        else:
            current_line = line
            #print(current_line)
            elements=line.split('\t')
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
            #print(current_spec)
            #current_spec = [read's UMI, strand (column2), chromosome(column3), LM_position(column4), cigar(column6), start_position=LM_position]
            if umi in known_umi_list:
                #print(umi)
                if bioinfo.flag_revcomp(strand) == 'forward':
                    #if strand is forward
                    #print(current_spec)
                    if re.match(r'^([0-9]|[1-9][0-9]|100)S.*', cigar):
                    #print(cigar)
                    #if S at begining of cigar
                        start_softclip=cigar.split('S')
                        #adjustS="##"
                        start_adjustS=int(start_softclip[0])
                        #print(adjustS)
                        start=start-start_adjustS
                        #print(start)
                        current_spec[5]=start
                        #print(current_spec)
                #print(current_spec)
                if bioinfo.flag_revcomp(strand) == 'reverse':
                #if strand is reverse
                    #print(current_spec)
                    if re.match(r'.*([0-9]|[1-9][0-9]|100)S$', cigar):
                    #if cigar from current_spec ends with "##s":
                        #print(cigar)
                        end_softclip=re.findall(r'.*([0-9]|[1-9][0-9]|100)S$', cigar)
                        #adjustS="##"
                        #print(end_softclip)
                        #print(end_softclip[0])
                        end_adjustS=int(end_softclip[0])
                        start=end_adjustS+start
                        current_spec[5]=start
                        #start position in current_spec=adjustS+LM_position
                #print(current_spec)
                    if re.findall(r'([0-9]|[1-9][0-9]|100)M.*', cigar):
                        #if "s' not in cigar from current_spec:
                        match=re.findall(r'([0-9]|[1-9][0-9]|100)M.*', cigar)
                        #print(match)
                        match_adjust=int(match[0])
                        #print(match_adjust)
                        start=match_adjust+start
                        current_spec[5]=start
                        #print(current_spec)
                        #start_position in current_spec = adjustM + LM_position
                    if re.findall(r'([0-9]|[1-9][0-9]|100)N.*', cigar):
                        #if cigar from current_spec has "N":
                        N_value=re.findall(r'([0-9]|[1-9][0-9]|100)N.*', cigar)
                        #print(N_value)
                        adjust_N=int(N_value[0])
                        start=adjust_N+start
                        current_spec[5]=start
                        #start_position in current_spec = adjustN + LM_position
                        #print(current_spec)
                    if re.findall(r'([0-9]|[1-9][0-9]|100)D.*', cigar):
                        D_value=re.findall(r'([0-9]|[1-9][0-9]|100)D.*', cigar)
                        adjust_D=int(D_value[0])
                        start=adjust_D+start
                        current_spec[5]=start
                        #print(current_spec)
                        #start_position in current_spec = adjustD + LM_position
                #print(current_spec)
                finaldict_empty=not bool(final_list)
                if finaldict_empty == True:
                    #print("first entry")
                    final_list=[current_spec[0], current_spec[1], current_spec[2], current_spec[5]]
                    #print(final_list)
                if finaldict_empty == False:
                    for item in final_list:
                        if current_spec[0] != final_list[0]:
                            final_list[].append(current_spec[0],current_spec[1],current_spec[2],current_spec[5])
                            #print("new umi")
                            #print(final_list)
                            current_spec=[]
                            break
                        else:
                            if current_spec[1] != final_list[1]:
                                final_list=[current_spec[0],current_spec[1],current_spec[2],current_spec[5]]
                                #print("new position")
                                #print(final_list)
                                current_spec=[]
                                break
                            else:
                                if current_spec[2] != final_list[2]:
                                    final_list=[current_spec[0],current_spec[1],current_spec[2],current_spec[5]]
                                    #print("new chrom")
                                    #print(final_list)
                                    current_spec=[]
                                    break
                                else:
                                    if current_spec[5] != final_list[3]:
                                        final_list=[current_spec[0],current_spec[1],current_spec[2],current_spec[5]]
                                        #print("new strand")
                                        #print(final_list)
                                        current_spec=[]
                                        break
                                    else:
                                        #print("copy")
                                        current_spec=[]
                                        break
    #print(final_list)

            #             print(final_list)
            #     * is UMI, start_position, chrom, strand in final array?
            #         * if yes:
            #             * empty current_line and current_spec
            #             * break
            #         * if no:
            #             * add info from current_spec to final array
            #             * write entire line from current_line to final sam file
            #             * empty current_line and current_spec
            #             * break
            # if UMI in current_spec array not in UMI_list:
            #   empty current_line and current_spec


            # current code assuming we writes adjusted position to final sam file