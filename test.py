import sys
import pysam
from pysam import VariantFile


def readfile(filename):
    bcf_in = VariantFile(filename, 'r')
    cnt = 1
    for rec in bcf_in.fetch():
        print(rec.pos)
        cnt += 1
        if cnt == 10:
            break


def main(argv):
    readfile(argv)


if __name__ == '__main__':
    main(sys.argv[1:])