from bam.bam_process import BAM
from vcf.vcf_process import VCF
from fasta.fasta_process import FASTA
import numpy as np


class DATAPROCESS:
    def __init__(self, vcf_filename, bam_filename, fasta_filename, data_filename):
        self.vcf_filename = vcf_filename
        self.bam_filename = bam_filename
        self.fasta_filename = fasta_filename
        self.data_filename = data_filename

    def dataproc(self):
        samples_data = []
        v = VCF()
        vcf_file = v.readfile(self.vcf_filename)
        ls_variant = v.varient_info(vcf_file)
        # print(ls_variant)
        b = BAM()
        bam_file = b.readfile(self.bam_filename)

        fa = FASTA()
        fasta_file = fa.readfile(self.fasta_filename)

        for rec in ls_variant:
            # 变异数据
            sample = rec[0].split('_')[0]
            chr = rec[0].split('_')[1]
            pos = rec[0].split('_')[2]
            ref = rec[1][0]
            ref = ref.lower()
            label = rec[2]
            variant_seq = b.pileup_column(bam_file, chr, int(pos), int(pos) + 1)
            if variant_seq is None:
                break

            variant_seq = [item.lower() for item in variant_seq]
            s_c_p = sample + '_' + chr + '_' + pos
            variant_sample = (s_c_p, (ref, tuple(variant_seq)), label)
            samples_data.append(variant_sample)

            # 非变异数据
            pos = int(pos) + 1
            while pos:
                normal_seq = b.pileup_column(bam_file, chr, pos, pos + 1)
                ref_base = fa.ref_atcg(fasta_file, chr, pos, pos + 1)
                ## 如果正常序列不存在就跳出
                if normal_seq == None or ref_base == None:
                    break

                normal_seq = [item.lower() for item in normal_seq]
                ref_base = ref_base.lower()
                norepeat = set(normal_seq)  ## 去重
                # print(norepeat)
                if len(norepeat) == 1 and list(norepeat)[0] == ref_base:
                    normal_s_c_p = sample + '_' + chr + '_' + str(pos)
                    normal_sample = (normal_s_c_p, (ref_base, tuple(normal_seq)), (0, 0))
                    # print(normal_sample)
                    samples_data.append(normal_sample)
                    break
                pos += 1

        np.save(self.data_filename, samples_data)
        # return samples_data


