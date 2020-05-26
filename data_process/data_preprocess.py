import pysam
from bam.bam_process import BAM
from vcf.vcf_process import VCF


class DATAPROCESS:
    def __init__(self, vcf_filename, bam_filename):
        self.vcf_filename = vcf_filename
        self.bam_filename = bam_filename

    def dataproc(self):
        sample_data = []
        v = VCF()
        vcf_file = v.readfile(self.vcf_filename)
        ls_variant = v.varient_info(vcf_file)
        print(ls_variant)
        b = BAM()
        bam_file = b.readfile(self.bam_filename)

        cnt = 1
        for rec in ls_variant:
            # print(rec)
            sample = rec[0].split('_')[0]
            chr = rec[0].split('_')[1]
            pos = rec[0].split('_')[2]
            ref = rec[1][0]
            label = rec[2]
            seq = b.pileup_column(bam_file, chr, int(pos), int(pos) + 1)
            s_c_p = sample + '_' + chr + '_' + pos
            sample = (s_c_p, (ref, tuple(seq)), label)
            sample_data.append(sample)
            print(sample)
            cnt += 1
            if cnt == 10:
                break




