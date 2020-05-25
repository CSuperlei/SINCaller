import pysam
from pysam import AlignmentFile

class BAM:
    def readfile(self, filename):
        bam_file = AlignmentFile(filename, 'rb')
        if None == bam_file:
            print("bam_file is empty")

        return bam_file

    def pileup_column(self, bam_file, chr_id, start, end):
        cnt = 1
        for rec in bam_file.pileup(chr_id, start, end):
            print(rec.pos)
            print(rec.get_mapping_qualities())
            print(rec.get_query_sequences(start))
            print(rec.get_query_positions())
            print(rec.reference_pos)
            cnt += 1
            if cnt == 10:
                break

