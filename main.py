import argparse
from vcf.vcf_process import VCF
from bam.bam_process import BAM


def args_func():
    parser = argparse.ArgumentParser(description="scSNV software")
    parser.add_argument('--vcf', '-v', help='vcf filename')
    parser.add_argument('--bam', '-b', help='bam filename')
    parser.add_argument('--fasta', '-fa', help='fasta, filename')
    args = parser.parse_args()
    return args


def main():
    args = args_func()

    # vcf_filename = args.vcf
    # v = VCF()
    # v.readfile(vcf_filename)

    bam_filename = args.bam
    b = BAM()
    bam_file = b.readfile(bam_filename)
    b.pileup_column(bam_file, 'chr1', 2956919, 2956922)


if __name__ == '__main__':
    main()



