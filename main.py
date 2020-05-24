import argparse
from vcf.vcf_process import VCF


def args_vcf():
    parser = argparse.ArgumentParser(description="VCF file")
    parser.add_argument('--vcf', '-v', help='vcf filename', required=True)
    args = parser.parse_args()
    return args


def main():
    args = args_vcf()
    # print(args.vcf)
    filename = args.vcf
    v = VCF(filename)
    v.readfile()


if __name__ == '__main__':
    main()



