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
                for tmp  in rec.pileups:
                    print(tmp.alignment.indel)
                return rec.get_query_sequences()

            elif rec.pos == end - 1:
                print('pos is not exist')
                return None
