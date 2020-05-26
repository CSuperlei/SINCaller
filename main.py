import argparse
from data_process.data_preprocess import DATAPROCESS


def args_func():
    parser = argparse.ArgumentParser(description="scSNV software")
    parser.add_argument('--vcf', '-v', help='vcf filename')
    parser.add_argument('--bam', '-b', help='bam filename')
    parser.add_argument('--fasta', '-fa', help='fasta, filename')
    parser.add_argument('--fastq', '-fq', help='fastq filename')
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


    vcf_filename = args.vcf
    bam_filename = args.bam
    d = DATAPROCESS(vcf_filename, bam_filename)
    samples_data = d.dataproc()
    print(samples_data)


if __name__ == '__main__':
    main()



