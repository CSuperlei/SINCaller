import numpy as np
from collections import Counter
import pysam
from pysam import AlignmentFile


class BAM:
    def readfile(self, filename):
        bam_file = AlignmentFile(filename, 'rb')
        if None == bam_file:
            print("bam_file is empty")

        return bam_file

    def pileup_column(self, bam_file, chr_id, start, end):
        for rec in bam_file.pileup(chr_id, start - 1, end - 1):  ## 索引从0开始
            if rec.pos == start - 1:
                # print(rec.get_mapping_qualities())
                # print(rec.get_query_sequences())
                # print(rec.get_query_positions())
                # print(rec.reference_pos)
                # print(dir(rec))
                # print(dir(rec.pileups))
                base_list = rec.get_query_sequences()
                indel_list = [int(tmp.indel) for tmp in rec.pileups]
                # print(indel_list)
                # for tmp in rec.pileups:
                #     print(dir(tmp))
                #     print(tmp.indel)
                #     print(dir(tmp.indel))
                #
                #     print(dir(tmp.alignment))
                #     print(tmp.alignment.reference_name)
                #     print(tmp.alignment.mapping_quality)
                def son_pilup(bam, chr, s, e):
                    for r in bam.pileup(chr, s, e):
                        if r.pos == s:
                            bl = r.get_query_sequences()
                            print(bl)
                            return bl

                re = '0'
                sum_indel_list = sum(indel_list)
                if sum_indel_list == 0:  ## 如果是SNP
                    bl = rec.get_query_sequences()
                elif sum_indel_list < 0:
                    bl = son_pilup(bam_file, chr_id, rec.pos + 1, rec.pos + 2)
                    indel_index = np.argmin(indel_list)
                    indel_value = np.min(indel_list)
                    re = self.fetch_row(bam_file, chr_id, rec.pos, rec.pos + 1, indel_index, indel_value)
                elif sum_indel_list > 0:
                    bl = son_pilup(bam_file, chr_id, rec.pos + 1, rec.pos + 2)
                    indel_value = np.max(indel_list)
                    indel_index = np.argmax(indel_list)
                    re = self.fetch_row(bam_file, chr_id, rec.pos, rec.pos + 1, indel_index, indel_value)

                # print(bl)
                base_ad = Counter(bl)
                ad = []  ## 计算不同等位基因数量
                for k, v in base_ad.items():
                    ad.append(v)
                dp = sum(ad) ## 总的映射深度
                ad = [str(i) for i in ad]
                ad = ",".join(ad) ## 每种等位基因的深度
                ad_dp = ad + '-' + str(dp)

                pileup_list = [base_list, indel_list, ad_dp, re]
                return pileup_list

            elif rec.pos == end - 1:
                print('pos is not exist')
                return None

    def fetch_row(self, bam_file, chr_id, start, end, index, indel_value):
        i = 0
        for rec in bam_file.fetch(chr_id, start - 1 , end - 1):
            if i != index:
                i += 1
                continue

            seq = list(rec.seq)
            reference = rec.get_reference_sequence()
            reference = list(reference)
            pairs = rec.get_aligned_pairs()
            if indel_value > 0: ## 插入
                indel_insertion = ""
                for item in pairs:
                    if start in item and None not in item:
                        ref = reference[item[0]]   ##找到indel插入的参考基因
                        for i in range(indel_value):
                            indel_insertion += seq[item[0] + i + 1]  ## 找到后边插入的基因是什么
                        indel_insertion = ref + indel_insertion
                        re = ref.upper() + '-' + indel_insertion.upper()
                        return re

            elif indel_value < 0:
                indel_deletion = ""
                for item in pairs:
                    if start in item and None not in item:
                        ref = reference[item[0]]  ## 找到indel缺失的参考基因
                        for i in range(-indel_value):
                            indel_deletion += reference[item[0] + i + 1]

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

    re = b.pileup_column(bam_file, 'chr1', 43785114, 43785115)
    print(re)

    re = b.pileup_column(bam_file, 'chr1', 44164156, 44164157)
    print(re)