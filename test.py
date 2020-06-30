# import pysam
# from pysam import VariantFile
# import argparse
# from bam.bam_process import BAM
# from fasta.fasta_process import FASTA
# import os
import argparse
from textwrap import dedent

# def readfile(filename):
#     bcf_in = VariantFile(filename, 'r')
#     cnt = 1
#     for rec in bcf_in.fetch():
#         for key, value in rec.samples.items():
#             print(key, value['AD'])
#         # print(rec.samples.iteritems())
#         print(rec.samples.get('AD'))
#         cnt += 1
#         if cnt == 10:
#             break


# def test_vcf():
#     parser = argparse.ArgumentParser(description="VCF file")
#     parser.add_argument('--vcf', '-v', help='vcf filename', required=True)
#     args = parser.parse_args()
#     return args


def main():
    # parser = argparse.ArgumentParser(description="scSNVIndel software")
    # parser.add_argument('--vcf', '-v', help='vcf filename')
    # parser.add_argument('--bam', '-b', help='bam filename')
    # parser.add_argument('--fasta', '-fa', help='fasta, filename')
    # parser.add_argument('--fastq', '-fq', help='fastq filename')
    # parser.add_argument('--gpus', '-g', help='gpu number')
    # parser.add_argument('--log', '-lo', help='log level')
    # parser.add_argument('--data', '-d', help='data filename')  ## 生成数据的名字
    # parser.add_argument('--data_model', '-dm',
    #                     help='data generator mode; mode 1 gernerates variant data, mode 2 generates normal data')
    # parser.add_argument('--load', '-ld', help='load filename')  ## 加载数据
    # parser.add_argument('--region', '-r', help='region test filename')  ## 加载测试区域数据
    # parser.add_argument('--test', '-tm',
    #                     help='test mode 1 is generator_test; mode 2 is batch test; mode 3 is random data')
    # parser.add_argument('--dc_origin', '-dco', help='data combine orgin')  ## 合并不同标签数据
    # parser.add_argument('--dc_target', '-dct', help='data combine target')  ## 生成不同标签数据
    # parser.add_argument('--mode', '-m',
    #                     help='mode 1 is training; mode 2 is tesing; mode 3 is generate data; mode 4 is combine data',
    #                     required=True)
    # args = parser.parse_args()
    # d = {'aa': 1, 'at': 2, 'ac': 3, 'ag': 4,
    #      'tt': 5, 'ta': 6, 'tc': 7, 'tg': 8,
    #      'cc': 9, 'ca': 10, 'ct': 11, 'cg': 12,
    #      'gg': 13, 'ga': 14, 'gc': 15, 'gt': 16}
    # docs = ['at at',
    #         'aa aa',
    #         'gt gt gt gt gt',
    #         'cc cc cc cc cc',
    #         'ag ag',
    #         'tt tt',
    #         'tc tc',
    #         'cc cc',
    #         'gc gc gc',
    #         'aa aa aa',
    #         'gg gg gt',
    #         'tt tt tt',
    #         'ag ag ag ag ag ag ag ag ag ag ag ag ag ag',
    #         'gg gg gg gg gg gg gg gg gg gg gg gg gg gg',
    #         'ga ga ga',
    #         'cc cc cc',
    #         ]
    #
    # for i, item in enumerate(docs):
    #     t = []
    #     for tmp in item.split(' '):
    #         r = d[tmp]
    #         t.append(r)
    #     docs[i] = t
    #     t = []
    #
    # print('docs')
    # print(docs)
    # labels = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
    # # labels = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    # # labels = [[0, 1, 2, 3], [4, 5, 6, 7]]
    # ohl = to_categorical(labels, num_classes=4)
    # print('ohl', ohl)
    # # integer encode the documents
    # vocab_size = 17
    # # encoded_docs = [one_hot(d, vocab_size) for d in docs]
    # # print('en', encoded_docs)
    # # pad documents to a max length of 4 words
    # max_length = 20
    # padded_docs = pad_sequences(docs, maxlen=max_length, padding='post')
    # print('padded_docs', padded_docs)
    # print('padded_docs shape', padded_docs.shape)

    # fasta_filename = args.fasta
    # bam_filename = args.bam
    # b = BAM()
    # bam_file = b.readfile(bam_filename)
    # print('hello')
    # fa = FASTA()
    # fasta_file = fa.readfile(fasta_filename)
    #
    # ref_base_indel = fa.ref_atcg(fasta_file, 'chr10', 2370295, 2370296)
    # print(ref_base_indel)
    # re = b.fetch_row(bam_file, 'chr10', 2370295, 2370296)
    # print(re)
    # re = b.pileup_column(bam_file, 'chr1', 10594829, 10594830)
    # print(re[0])
    # print(re[1])
    output_file = open('./out.txt', 'w')
    def output(string):
        print(string, file=output_file)

    output(dedent(
            """\
            ##fileformat=VCFv4.2
            ##FILTER=<ID=LowQual,Description="Low quality">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
            ##scSNVIndelCommandLine=<ID=HaplotypeCaller,CommandLine="python main.py  -ld /home/cailei/bio_project/nbCNV/predict_result/SRR053608_predict_data_chr1.npy -ov /home/cailei/bio_project/nbCNV/scSNVIndel_vcf/SRR053608_predict_chr1.vcf -g 56 -lo 3 -m 5">
            ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
            ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
            ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
            ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
            ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
            ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
            ##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
            ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
            ##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
            ##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
            ##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
            ##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
            ##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
            ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
            ##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
            ##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">"""
        )

    )

    output(dedent('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSRR05107'))
    output('a\tb\tc\td\te\tf\t')
    output_file.close()

if __name__ == '__main__':
    main()