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
            s_c_p = sample + '_' + rec.chrom + '_' + str(rec.pos)
            self.ls.append((s_c_p, (rec.ref, rec.alts[0]), label))
            print(self.ls)
            cnt += 1
            if cnt == 10:
                break








