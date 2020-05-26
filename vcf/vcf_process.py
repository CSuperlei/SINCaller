import pysam
from pysam import VariantFile

class VCF:
    def readfile(self, filename):
        bcf_in = VariantFile(filename, 'r')
        if None == bcf_in:
            print('vcf file is empty')

        return bcf_in


    def varient_info(self, vcf_file):
        ls = []
        for rec in vcf_file.fetch():
            sample = ''
            label = ()
            for key, value in rec.samples.items():
                sample = key
                label = value['GT']
                break
            s_c_p = sample + '_' + rec.chrom + '_' + str(rec.pos)
            ls.append((s_c_p, (rec.ref, rec.alts[0]), label))

        return ls








