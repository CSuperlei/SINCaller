import sys
from vcf.vcf_process import VariantFile

def main(argv):
    v = VariantFile(argv)
    v.readfile()

if __name__ == '__main__':
    main(sys.argv[1:])



