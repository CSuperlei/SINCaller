import sys
from vcf.vcf_process import VCF

def main(argv):
    v = VCF(argv)
    v.readfile()

if __name__ == '__main__':
    main(sys.argv[1:])



