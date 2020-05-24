import pysam
from pysam import VariantFile

class VCF:
    def __init__(self, filename):
        self.filename = filename
        self.ls = []

    def readfile(self):
        bcf_in = VariantFile(self.filename, 'r')
        cnt = 1
        for rec in bcf_in.fetch():
            sample = ''
            label = ()
            for key, value in rec.samples.items():
                sample = key
                label = value['GT']
                break
            self.ls.append((sample, rec.chrom, rec.pos, (rec.ref, rec.alts), label))
            print(self.ls)
            cnt += 1
            if cnt == 10:
                break








