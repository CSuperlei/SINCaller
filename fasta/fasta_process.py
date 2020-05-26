import pysam
from pysam import FastaFile


class FASTA:
    def readfile(self, filename):
        fasta_file = FastaFile(filename)
        if None == fasta_file:
            print("fasta_file is empty")

        return fasta_file

    def ref_atcg(self, fasta_file, chr_id, start, end):
        for rec in fasta_file.fetch(chr_id, start - 1, end - 1):
            # print(rec)  ## 返回某个位点参考基因
            return rec
            break
