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
                print(indel_list)
                # for tmp in rec.pileups:
                #     print(dir(tmp))
                #     print(tmp.indel)
                #     print(dir(tmp.indel))
                #
                #     print(dir(tmp.alignment))
                #     print(tmp.alignment.reference_name)
                #     print(tmp.alignment.mapping_quality)

                pileup_list = [base_list, indel_list]
                return pileup_list

            elif rec.pos == end - 1:
                print('pos is not exist')
                return None


    def fetch_row(self, bam_file, chr_id, start, end):
        for rec in bam_file.fetch(chr_id, start - 1 , end - 1):
            print(type(rec.get_reference_positions))
            offset = int(start - 1) - int(rec.get_reference_positions[0])
            seq = list(rec.seq)
            print(seq)
            print(seq[offset])

