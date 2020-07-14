import numpy as np
from collections import Counter
import pysam
from pysam import AlignmentFile
from fasta.fasta_process import FASTA


class BAM:
    def readfile(self, filename):
        bam_file = AlignmentFile(filename, 'rb')
        if None == bam_file:
            print("bam_file is empty")

        return bam_file

    def pileup_column(self, bam_file, chr_id, start, end, fasta_file):
        for rec in bam_file.pileup(chr_id, start - 1, end - 1, stepper='all', ignore_overlaps=True):  ## 索引从0开始
            if rec.pos == start - 1:
                base_list = rec.get_query_sequences()
                indel_list = [int(tmp.indel) for tmp in rec.pileups]

                sum_indel_list = sum(indel_list)
                if sum_indel_list == 0:  ## 如果是SNP
                    re = '0'
                elif sum_indel_list < 0:
                    indel_index = np.argmin(indel_list)
                    indel_value = np.min(indel_list)
                    re = self.fetch_row(bam_file, chr_id, rec.pos, rec.pos + 1, indel_index, indel_value, fasta_file)
                elif sum_indel_list > 0:
                    indel_value = np.max(indel_list)
                    indel_index = np.argmax(indel_list)
                    re = self.fetch_row(bam_file, chr_id, rec.pos, rec.pos + 1, indel_index, indel_value, fasta_file)

                base_ad = Counter(indel_list)
                ad = []  ## 计算不同等位基因数量
                for k, v in base_ad.items():
                    if k != 0:
                        ad.append(v)
                dp = len(indel_list) ## 总的映射深度
                d = dp - sum(ad)  ## 与参考基因相同的数量
                if sum(indel_list) != 0:
                    ad = [str(i) for i in ad]
                    ad = ",".join(ad) ## 每种等位基因的深度
                    ad = str(d) + ',' + ad
                    ad_dp = ad + '-' + str(dp)
                else:
                    ad_dp = str(dp) + '-' + str(dp)

                pileup_list = [base_list, indel_list, ad_dp, re]
                return pileup_list

            elif rec.pos == end - 1:
                print('pos is not exist')
                return None


    def fetch_row(self, bam_file, chr_id, start, end, index, indel_value, fasta_file):
        fa = FASTA()
        i = 0
        for rec in bam_file.fetch(chr_id, start - 1 , end - 1, multiple_iterators=True, until_eof=True):
            if i != index:
                i += 1  ## 找发生缺失的那条read
                continue

            seq = list(rec.seq)
            reference = rec.get_reference_sequence()
            reference = list(reference)
            pairs = rec.get_aligned_pairs()
            if indel_value > 0: ## 插入
                indel_insertion = ""
                for item in pairs:
                    if start in item and None not in item:
                        # ref = reference[item[0]]   ##找到indel插入的参考基因
                        ref = fa.ref_atcg(fasta_file, item[1], item[1] + 1)
                        for i in range(indel_value):
                            indel_insertion += seq[item[0] + i + 1]  ## 找到后边插入的基因是什么
                        indel_insertion = ref + indel_insertion
                        re = ref.upper() + '-' + indel_insertion.upper()
                        return re

            elif indel_value < 0:
                indel_deletion = ""
                for item in pairs:
                    if start in item and None not in item:
                        # ref = reference[item[0]]  ## 找到indel缺失的参考基因
                        ref = fa.ref_atcg(fasta_file, item[1], item[1] + 1)
                        for i in range(-indel_value):
                            # indel_deletion += reference[item[0] + i + 1] ## 找到缺失的参考基因是什么
                            indel_deletion += fa.ref_atcg(fasta_file, item[1] + i + 1, item[1] + i + 2);

                        indel_deletion = ref + indel_deletion
                        re = indel_deletion.upper() + '-' + ref.upper()
                        return re


if __name__ == '__main__':
    f_d = '/home/cailei/bio_project/nbCNV/bam/SRR052047/SRR052047_rmduplicate.bam'
    b = BAM()
    bam_file = b.readfile(f_d)
    # re = b.fetch_row(bam_file, 'chr1', 31988447, 31988448, 2, -1)
    # print(re)
    #
    # re = b.fetch_row(bam_file, 'chr1', 43785114, 43785115, 1, 1)
    # print(re)

    # re = b.pileup_column(bam_file, 'chr1', 96689497, 96689498)
    # print(re)
    #
    # re = b.pileup_column(bam_file, 'chr1', 44164156, 44164157)
    # print(re)

    column = b.pileup_column(bam_file, 'chr1', 45571770, 45571771)
    print(column)

    column = b.pileup_column(bam_file, 'chr1', 45571771, 45571772)
    print(column)

    column = b.pileup_column(bam_file, 'chr1', 45571772, 45571773)
    print(column)
