Write up a strategy for writing a Reference Based PCR Duplicate Removal tool. That is, given a sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory. You should not write any code for this portion of the assignment. Be sure to:

- Define the problem:
    * We need to identify and remove PCR duplicates created during the laboratory process so we retain only the original copies. The input will be a sam file of mapped reads and the output will be a file with only a single copy of each read. 
- Write examples:
    - Include a properly formatted input sam file (an example sam file is included in the repo for your reference)
    - Include a properly formatted expected output sam file
    - Cover several different cases of things that are and are not PCR duplicates
    - It may be helpful to create a "unit test folder", that contains it's own readme.md describing the test cases in your unit test samfile

- Develop your algorithm using pseudocode:

* create UMI_list from provided UMI file
* create empty array (final_list): will hold UMI, position, chrom, strand for each read and will compare future reads for duplicates
* create temporary array (current_spec): will hold current UMI, position, chrom, strand
* with open sam file:
    * if line starts with "@" (sam file header):
        * write line to final sam file
    * if line doesn't start with "@":
        * if line is empty "":
            * break
        * else:
            * current_line = entire read
            * UMI = column1, awk statement to find 8th item separated by :
            * current_spec = [read's UMI, strand (column2), chromosome(column3), LM_position(column4), cigar(column6), start_position=LM_position]
            * if UMI in current_spec array is in UMI_list:
                * if strand_function applied to strand in current_spec = mapped:
                    * if bitwise_revcomp = false:
                        * if cigar from current_spec starts with "##s":
                            * adjustS="##"
                            * start_position in current_spec=LM_position-adjustS
                        * if cigar from current_spec ends with "##s":
                            * ?
                        * if "s' not in cigar from current_spec:
                            * continue
                    * if bitwise_revcomp applied to strand in current_spec = true:
                        * if cigar from current_spec starts with "##s":
                            * adjustS="##"
                            * start_position in current_spec=adjustS+LM_position
                        * if cigar from current_spec ends with "##s":
                            * continue
                        * if "s' not in cigar from current_spec:
                            * adjustM= ## before m
                            * start_position in current_spec = adjustM + LM_position
                * is UMI, start_position, chrom, strand in final array?
                    * if yes:
                        * empty current_line and current_spec
                        * break
                    * if no:
                        * add info from current_spec to final array
                        * write entire line from current_line to final sam file
                        * empty current_line and current_spec
                        * break
                * if strand_function == unmapped:
                    * empty current_line and current_spec
                    * break
            * if UMI in current_spec array not in UMI_list:
                * empty current_line and current_spec
                * break

- Determine high level functions

def strand_function:
    "determine if bitwise flag in sam file indicates read is mapped (true) or unmapped (false)"
    if ((flag & 4) == 4):
        mapped = false
        return 'unmapped'
    if ((flag & 4) != 4):
        mapped = true
        return 'mapped'
def bitwise_revcomp:
    "determine if bitwise flag in sam file indicates read is reverse complement"
    if ((flag & 16) == 16): 
    rev_comp = true
    return 'true'

assert strand_function("")
assert bitwise_revcom("")



For this portion of the assignment, you should design your algorithm for single-end data, with 96 known UMIs. A list of UMIs can also be found in the repo. UMI information will be in the QNAME, like so: NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT. Discard any UMIs with errors (or error correct, if you're feeling ambitious).