import sys
import pysam
from pysam import VariantFile
import argparse

def readfile(filename):
    bcf_in = VariantFile(filename, 'r')
    cnt = 1
    for rec in bcf_in.fetch():
        # for key, value in rec.samples.iteritems():
        #     print(key, value['GT'])
        print(rec.samples['GT'])
        cnt += 1
        if cnt == 10:
            break


def test_vcf():
    parser = argparse.ArgumentParser(description="VCF file")
    parser.add_argument('--vcf', '-v', help='vcf filename', required=True)
    args = parser.parse_args()
    return args


def main():
    args = test_vcf()
    print(args.vcf)
    filename = args.vcf
    readfile(filename)


if __name__ == '__main__':
    main()