import numpy as np
from collections import Counter
from vcf.vcf_process import VCF
from bam.bam_process import BAM


class GVCF:
    def __init__(self, load_filename, vcf_filename, fastai_name):
        self.load_filename = load_filename ## 加载模型预测数据
        self.vcf_filename = vcf_filename   ## 生成vcf文件
        self.fastai_name = fastai_name     ## 加载染色体长度文件
        self.diction = {
            1: 'aa', 2: 'ac', 3: 'ag', 4: 'at', 5: 'ad',
            6: 'cc', 7: 'ca', 8: 'cg', 9: 'ct', 10: 'cd',
            11: 'gg', 12: 'ga', 13:'gc', 14:'gt', 15: 'gd',
            16: 'tt', 17: 'ta', 18: 'tc', 19: 'tg', 20: 'td'
        }

    def __int_to_str(self, i):
        s = self.diction[i]
        return s

    def __filter(self, p):
        if p < 0.5:
            return 'LowQual'
        else:
            return 'PASS'

    def __gt(self, g):
        if g == 0:
            return ''
        if g == 1:
            return '0/1'
        if g == 2:
            return '1/1'

    def generation_vcf(self):
        v = VCF()
        title = """\
                ##fileformat=VCFv4.2
                ##FILTER=<ID=LowQual,Description="Low quality">
                ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
                ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
                ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
                ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                ##FORMAT=<ID=PL,Number=G,Type=Float,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
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

        data = np.load(self.load_filename, allow_pickle=True)
        sample_name = data[0][0].split('_')[0]
        ### title, fastai_name, sample_name, vcf_filename_file
        v.gernate_vcf_title(title, self.fastai_name, sample_name, self.vcf_filename)

        i = 0
        while i < len(data):
            SAMPLE = data[i][0].split('_')[0]
            CHROM = data[i][0].split('_')[1]
            POS = data[i][0].split('_')[2]

            indel_value = data[i][0].split('_')[3]
            AD = data[i][0].split('_')[4].split('-')[0]
            DP = data[i][0].split('_')[4].split('-')[1]
            if int(DP) < 2:
                i += 1
                continue
            REF_ALT = data[i][0].split('_')[5]
            ## 预处理
            ## 判断是否发生了碱基变化
            base_pair = data[i][1][0]
            base_pair_prob = data[i][1][1]

            ## 发生了indel变异
            indel_pair = data[i][2][0]
            indel_pair_pro = data[i][2][1]

            genotype_pair = data[i][3][0]
            genotype_pair_pro = data[i][3][1]

            PL = str(round(1 - data[i][4][0], 3)) + ',' + str(round(1 - data[i][4][1], 3)) + ',' + str(round(1 - data[i][4][2],3))

            REF = self.__int_to_str(base_pair)[0]
            ALT = self.__int_to_str(base_pair)[1]

            sum_e = np.exp(base_pair_prob) + np.exp(indel_pair_pro) + np.exp(genotype_pair_pro)

            GT = self.__gt(genotype_pair)

            ## SNP 变异
            if REF != ALT:
                ID = '.'
                REF = REF.upper()
                ALT = ALT.upper()
                QUAL = round(-10 * np.log10(1 - (np.exp(base_pair_prob) + np.exp(genotype_pair_pro)) / sum_e), 3)
                if QUAL != 4.77:
                    i += 1
                    continue
                FILTER = self.__filter(QUAL)
                INFO = '.'
                FORMAT = 'GT:AD:DP:GQ:PL'
                VALUE = GT + ':' + AD + ':' + DP + ':' + str(QUAL) + ':' + PL
                v.generate_vcf_content(self.vcf_filename, CHROM, POS, ID, REF, ALT, str(QUAL), FILTER, INFO, FORMAT, VALUE)
                i += 1
            ## Indel 变异
            elif REF == ALT and indel_pair != 0:
                ID = '.'
                if REF_ALT == 'None' or REF_ALT == '0' or REF_ALT == None:
                    i = i + 1
                    continue
                REF = REF_ALT.split('-')[0]
                ALT = REF_ALT.split('-')[1]
                if QUAL <= 4.77:
                    i += 1
                    continue
                QUAL = round(-10 * np.log10(1 - (np.exp(indel_pair_pro) + np.exp(genotype_pair_pro)) / sum_e), 3)
                FILTER = self.__filter(QUAL)
                INFO = '.'
                FORMAT = 'GT:AD:DP:GQ:PL'
                VALUE = GT + ':' + AD + ':' + DP + ':' + str(QUAL) + ':' + PL
                v.generate_vcf_content(self.vcf_filename, CHROM, POS, ID, REF, ALT, str(QUAL), FILTER, INFO, FORMAT, VALUE)
                i += int(indel_value)
            else:
                i += 1

