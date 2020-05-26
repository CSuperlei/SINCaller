import pysam
from pysam import FastxFile


class FASTQ:
    def readfile(self, filename):
        fastq_file = FastxFile(filename)   ## 必须是.gz格式
        if None == fastq_file:
            print("fastq_file is empty")

        return fastq_file

    def seq_atcg(self, fastq_file):
        cnt = 1
        for rec in fastq_file:
            print(rec.name)
            print(rec.sequence)
            print(rec.comment)
            print(rec.quality)
            cnt += 1
            if cnt == 10:
                break
