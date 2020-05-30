import os
import argparse
import numpy as np
from data_process.data_preprocess import DATAPROCESS
from model.train_process import training

def args_func():
    parser = argparse.ArgumentParser(description="scSNV software")
    parser.add_argument('--vcf', '-v', help='vcf filename')
    parser.add_argument('--bam', '-b', help='bam filename')
    parser.add_argument('--fasta', '-fa', help='fasta, filename')
    parser.add_argument('--fastq', '-fq', help='fastq filename')
    parser.add_argument('--gpus', '-g', help='gpu number')
    args = parser.parse_args()
    return args


def main():
    args = args_func()

    # vcf_filename = args.vcf
    # v = VCF()
    # v.readfile(vcf_filename)

    # bam_filename = args.bam
    #     # b = BAM()
    #     # bam_file = b.readfile(bam_filename)
    #     # b.pileup_column(bam_file, 'chr1', 2956920, 2956921)

    # fasta_filename = args.fasta
    # a = FASTA()
    # fasta_file = a.readfile(fasta_filename)
    # a.ref_atcg(fasta_file, 'chr1', 2956920, 2957922)

    # fastq_filename = args.fastq
    # q = FASTQ()
    # fastq_file = q.readfile(fastq_filename)
    # q.seq_atcg(fastq_file)

    gpus = args.gpus
    gpus = ",".join(list(str(gpus)))
    os.environ["CUDA_VISIBLE_DEVICES"] = gpus

    vcf_filename = args.vcf
    bam_filename = args.bam
    fasta_filename = args.fasta
    d = DATAPROCESS(vcf_filename, bam_filename, fasta_filename)
    samples_data = d.dataproc()
    print(samples_data)
    sep = int(np.floor(len(samples_data) * 0.8))
    samples_train_data = samples_data[:sep]
    samples_test_data = samples_data[sep:]

    sendin = training(samples_train_data, samples_test_data)
    print(sendin)


if __name__ == '__main__':
    main()



