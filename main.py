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
    parser.add_argument('--log', '-lo', help='log level')
    parser.add_argument('--data', '-d', help='data filename')
    parser.add_argument('--load', '-ld', help='load filename')
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
    print('gpus', gpus)
    os.environ["CUDA_VISIBLE_DEVICES"] = gpus

    log_level = args.log
    print('log_level', log_level)
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = log_level

    vcf_filename = args.vcf
    bam_filename = args.bam
    fasta_filename = args.fasta
    data_filename = args.data_filename
    load_filename = args.load_filename

    if load_filename is not None:
        samples_data = np.load(load_filename)
    elif data_filename is not None:
        d = DATAPROCESS(vcf_filename, bam_filename, fasta_filename, data_filename)
        samples_data = d.dataproc()

    if samples_data is not None:
        # print(samples_data)
        sep = int(np.floor(len(samples_data) * 0.8))
        samples_train_data = samples_data[:sep]
        samples_test_data = samples_data[sep:]
        sendin = training(samples_train_data, samples_test_data)
        print(sendin)
    else:
        print('data is empty')


if __name__ == '__main__':
    main()



