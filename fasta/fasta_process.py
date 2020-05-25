import pysam
from pysam import FastaFile


class FASTA:
    def readfile(self, filename):
        fasta_file = FastaFile(filename, 'rb')
        if None == fasta_file:
            print("fasta_file is empty")

        return fasta_file

    def ref_atcg(self, fasta_file, start, end):
        for rec in fasta_file.fetch(start=start, end=end):
            print(rec)
