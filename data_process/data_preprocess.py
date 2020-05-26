from bam.bam_process import BAM
from vcf.vcf_process import VCF
from fasta.fasta_process import FASTA


class DATAPROCESS:
    def __init__(self, vcf_filename, bam_filename, fasta_filename):
        self.vcf_filename = vcf_filename
        self.bam_filename = bam_filename
        self.fasta_filename = fasta_filename

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

        cnt = 1
        for rec in ls_variant:
            # 变异数据
            sample = rec[0].split('_')[0]
            chr = rec[0].split('_')[1]
            pos = rec[0].split('_')[2]
            ref = rec[1][0]
            label = rec[2]
            variant_seq = b.pileup_column(bam_file, chr, int(pos), int(pos) + 1)
            s_c_p = sample + '_' + chr + '_' + pos
            variant_sample = (s_c_p, (ref, tuple(variant_seq)), label)
            samples_data.append(variant_sample)

            # 非变异数据
            pos = int(pos) + 1
            while pos:
                normal_seq = b.pileup_column(bam_file, chr, pos, pos + 1)
                ref_base = fa.ref_atcg(fasta_file, chr, pos, pos + 1)
                norepeat = set(normal_seq)
                if len(norepeat) == 1 and norepeat[0] == ref_base:
                    normal_s_c_p = sample + '_' + chr + '_' + str(pos)
                    normal_sample = (normal_s_c_p, (ref_base, tuple(normal_seq)), (0, 0))
                    samples_data.append(normal_sample)
                    break
                pos += 1

            cnt += 1
            if cnt == 10:
                break

        return samples_data


