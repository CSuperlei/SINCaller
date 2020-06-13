from bam.bam_process import BAM
from vcf.vcf_process import VCF
from fasta.fasta_process import FASTA
from data_process.region_process import TREGION
import numpy as np


class DATAPROCESS:
    def __init__(self, vcf_filename, bam_filename=None, fasta_filename=None, data_filename=None, region_filename=None):
        self.vcf_filename = vcf_filename
        self.bam_filename = bam_filename
        self.fasta_filename = fasta_filename
        self.data_filename = data_filename
        self.region_filename = region_filename
        self.l = {'a': 26, 'c': 27, 'g': 28, 't': 29, 'd': 30,
                  'aa': 31, 'ac': 32, 'ag': 33, 'at': 34, 'ad': 35,
                  'cc': 36, 'ca': 37, 'cg': 38, 'ct': 39, 'cd': 40,
                  'gg': 41, 'ga': 42, 'gc': 43, 'gt': 44, 'gd': 45,
                  'tt': 46, 'ta': 47, 'tc': 48, 'tg': 49, 'td': 50,
                  }

    def __str_to_int(self, s):
        r = self.l[s]
        return r

    ## 生成已知label数据
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

            variant_seq = ['d' if item == '' else item for item in variant_seq]
            variant_seq = [item.lower() for item in variant_seq]
            s_c_p = sample + '_' + chr + '_' + pos
            variant_seq = [ref + i for i in variant_seq]
            variant_seq = [self.__str_to_int(i) for i in variant_seq]
            # variant_sample = (s_c_p, (ref, tuple(variant_seq)), label)
            variant_sample = (s_c_p, variant_seq, 1)
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
                    normal_seq = [ref_base + i for i in normal_seq]
                    normal_seq = [self.__str_to_int(i) for i in normal_seq]
                    # normal_sample = (normal_s_c_p, (ref_base, tuple(normal_seq)), (0, 0))
                    normal_sample = (normal_s_c_p, normal_seq, 0)
                    # print(normal_sample)
                    samples_data.append(normal_sample)
                    break
                pos += 1

        np.save(self.data_filename, samples_data)
        return samples_data


    ## 生成未知label数据
    def test_pos(self):
        samples_data = []

        b = BAM()
        bam_file = b.readfile(self.bam_filename)

        fa = FASTA()
        fasta_file = fa.readfile(self.fasta_filename)

        t = TREGION()
        region_file = t.readfile(self.region_filename)
        region = t.region_info(region_file)

        for rec in region:
            sample = rec[0]
            chr = rec[1]
            l_pos = rec[2]
            r_pos = rec[3]
            for pos in range(int(l_pos), int(r_pos) + 1):
                seq = b.pileup_column(bam_file, chr, pos, pos + 1)
                ref_base = fa.ref_atcg(fasta_file, chr, pos, pos + 1)

                if ref_base is None or seq is None:
                    break

                seq = ['d' if item == '' else item for item in seq]
                seq = [item.lower() for item in seq]
                ref_base = ref_base.lower()

                s_c_p = sample + '_' + chr + '_' + str(pos)
                samples_data.append((s_c_p, (ref_base, tuple(seq))))

        np.save(self.data_filename, samples_data)
        return samples_data




