#!/bin/env python3
from optparse import OptionParser
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor

def main():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", "--ref", action="store", dest="ref", help="reference (fasta)", default="")
    parser.add_option("-s", "--snps", action="store", dest="snps", help="snp table (CSV)", default="")
    parser.add_option("-v", "--vcf", action="store", dest="vcf", help="vcf file", default="")
    parser.add_option("-i", "--id", action="store", dest="sample", help="Sample ID", default="")
    parser.add_option("-o", "--output", action="store", dest="output", help="output file", default="Output.fasta")
    return parser.parse_args()

def find_seq(snplist, loc, alt, ref):
    seq = []
    for i in range(0,len(snplist)):
        if snplist[i] in loc:
            j = loc.index(snplist[i])
            seq.append(alt[j])
        else:
            seq.append(ref.seq[int(snplist[i])-1])
    return seq

if __name__ == "__main__":
    options, args = main()
    output = open(options.output, 'w')
    idd = str(options.sample)

    # Read reference sequence
    ref = SeqIO.read(options.ref, "fasta")

    # Read in SNP info
    f = open(options.vcf,"r")
    loc = []
    alt = []
    for line in f:
        if "#" not in line:
            fields = line.rstrip().split("\t")
            loc.append(fields[1])
            alt.append(fields[4])
    f.close()

    # Read in SNP locations
    p = open(options.snps,"r")
    snplist = []
    for line in p:
        snplist.append(line.rstrip())
    p.close()

    with ThreadPoolExecutor() as executor:
        futures = []
        chunk_size = 1000
        for i in range(0, len(snplist), chunk_size):
            chunk = snplist[i:i+chunk_size]
            futures.append(executor.submit(find_seq, chunk, loc, alt, ref))
        seq_chunks = [future.result() for future in futures]

    seq = [item for sublist in seq_chunks for item in sublist]
    record = '>'+idd+'\n'+(''.join(seq))
    output.write(record)
    output.close()
