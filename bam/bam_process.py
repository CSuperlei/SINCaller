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
                base_list = rec.get_query_sequences()
                # print(rec.get_query_positions())
                # print(rec.reference_pos)
                # print(dir(rec))
                # print(dir(rec.pileups))
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
                sum_indel_list = sum(indel_list)
                if sum_indel_list == 0:
                    base_ad = Counter(base_list)
                    ad = []  ## 计算不同等位基因数量
                    for k, v in base_ad.items():
                        if k != '':
                            ad.append(v)
                    ad = ",".join(ad)
                    dp = rec.n ## 计算Coverage深度

                if sum_indel_list < 0:
                    deletion = np.min(indel_list)

                if sum_indel_list > 0:
                    insertion = np.max(indel_list)
                    index = np.argmax(indel_list)
                    re = self.fetch_row(bam_file, chr_id, rec.pos, rec.pos + insertion)
                    insert_seq = re[index]


                pileup_list = [base_list, indel_list, ad, dp]
                return pileup_list

            elif rec.pos == end - 1:
                print('pos is not exist')
                return None

    def fetch_row(self, bam_file, chr_id, start, end):
        re = []
        for rec in bam_file.fetch(chr_id, start - 1 , end - 1):
            print(start - 1)
            print(end - 1)
            print(list(rec.get_reference_positions()))
            print((rec.get_reference_sequence()))
            ## 求出当前位点到序列起始位点的长度
            offset = int(start - 1) - int(rec.get_reference_positions()[0])
            seq = list(rec.seq)
            print(rec.seq)
            print(rec.get_aligned_pairs())
            # print(seq[offset])
            ## 求得匹配到该点的序列
            re.append(seq[offset])

        return re

if __name__ == '__main__':
    f_d = '/home/cailei/bio_project/nbCNV/bam/SRR052047/SRR052047_rmduplicate.bam'
    b = BAM()
    bam_file = b.readfile(f_d)
    re = b.fetch_row(bam_file, 'chr1', 31988412, 31988476)
    print(re)

    re = b.fetch_row(bam_file, 'chr1', 43785086, 43785157)
    print(re)

