import pysam
from pysam import VariantFile

class VCF:
    def __init__(self, filename):
        self.filename = filename
        self.ls = ()

    def readfile(self):
        bcf_in = VariantFile(self.filename, 'r')
        for rec in bcf_in.fetch():
            self.ls.append(rec.chrom, rec.pos, (rec.ref, rec.alt))
            print(self.ls)








