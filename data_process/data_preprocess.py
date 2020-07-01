from bam.bam_process import BAM
from vcf.vcf_process import VCF
from fasta.fasta_process import FASTA
from data_process.region_process import TREGION
import numpy as np
from scipy import stats
from collections import Counter


class DATAPROCESS:
    def __init__(self, mode=1, vcf_filename=None, bam_filename=None, fasta_filename=None, data_filename=None, region_filename=None, padded_maxlen=78):
        self.vcf_filename = vcf_filename
        self.bam_filename = bam_filename
        self.fasta_filename = fasta_filename
        self.data_filename = data_filename
        self.region_filename = region_filename
        self.padded_maxlen = padded_maxlen
        self.mode = mode
        self.l = {
                  'aa': 1, 'ac': 2, 'ag': 3, 'at': 4, 'ad': 5,
                  'cc': 6, 'ca': 7, 'cg': 8, 'ct': 9, 'cd': 10,
                  'gg': 11, 'ga': 12, 'gc': 13, 'gt': 14, 'gd': 15,
                  'tt': 16, 'ta': 17, 'tc': 18, 'tg': 19, 'td': 20,
                   'a': 21, 'a+': 22, 'a-': 23,
                   'c': 24, 'c+': 25, 'c-': 26,
                   'g': 27, 'g+': 28, 'g-': 29,
                   't': 30, 't+': 31, 't-': 32
                  }

    def __str_to_int(self, s):
        r = self.l[s]
        return r

    def __indel_to_int(self, item):
        if item < 0:
            return 33
        elif item == 0:
            return 34
        elif item > 0:
            return 35

    def __padded_fill(self, data=None, padded_len=78):
        data_len = len(data)
        fill_zero = padded_len - data_len
        zero_list = [0 for i in range(fill_zero)]
        padded_result = data + zero_list
        return padded_result

    ## 生成已知label数据
    def dataproc(self):
        samples_data = []
        v = VCF()
        vcf_file = v.readfile(self.vcf_filename)
        ls_variant = v.varient_info(vcf_file)
        b = BAM()
        bam_file = b.readfile(self.bam_filename)

        fa = FASTA()
        fasta_file = fa.readfile(self.fasta_filename)

        for rec in ls_variant:
            # 变异数据
            sample = rec[0].split('_')[0]
            chr = rec[0].split('_')[1]
            pos = rec[0].split('_')[2]
            ref_alts = rec[1][0]
            alts = rec[1][1]
            if pos.isdigit() == False:
                continue
            pos = int(pos)
            if int(self.mode) == 1: ## 生成变异数据
                ref_base = fa.ref_atcg(fasta_file, chr, pos, pos + 1)
                ref = ref_base.lower()
                label = rec[2]
                pileup_list = b.pileup_column(bam_file, chr, pos, pos + 1)
                if pileup_list is None or pileup_list[0] is None or pileup_list[1] is None or ref_base is None or alts == '*':
                    continue

                ## 处理该位点得碱基序列
                variant_seq = ['d' if item == '' else item for item in pileup_list[0]]
                lower_list = [item.lower() for item in variant_seq]
                ref_var_list = [ref + i for i in lower_list]
                ref_var_list = [self.__str_to_int(i) for i in ref_var_list]

                ## 处理该位点的indel序列
                indel_list = pileup_list[1]
                indel_seq = [self.__indel_to_int(i) for i in indel_list]

                ## 处理genotype
                indel_sum = sum(indel_list)
                ref_base_genotype = fa.ref_atcg(fasta_file, chr, pos, pos + 1)  ## 读下一个位置
                ref_base_genotype = ref_base_genotype.lower()
                ## snp基因型处理
                if indel_sum == 0:
                    genotype_list = ref_var_list

                    ## 计算碱基改变系数
                    # base_coefficient = [1 if ref != i else 0 for i in lower_list]
                    # base_coefficient = round(sum(base_coefficient) / len(base_coefficient), 3)

                ## indel 基因型处理
                elif indel_sum > 0:
                    ## 插入变异
                    genotype_list = [ref_base_genotype + '+' if i > 0 else ref_base_genotype for i in indel_list]
                    genotype_list = [self.__str_to_int(i) for i in genotype_list]

                    ## 碱基改变系数
                    # base_coefficient = [1 if i > 0 else 0 for i in indel_list]
                    # base_coefficient = round(sum(base_coefficient) / len(indel_list), 3)

                elif indel_sum < 0:
                    '''
                    在indel缺失有一种特殊情况,就是发生缺失又发生碱基改变，在产生indel序列时
                    出现[0,0,0,-1]其实0的后一个位置与参考基因不同,只是碱基发生了改变，所以会出现0的情况，
                    本质上还是发生了indel缺失变异
                    '''
                    ## 如果是缺失indel要特判一下
                    ref_base_indel_next = fa.ref_atcg(fasta_file, chr, pos + 1, pos + 2) ## 取indel缺失位置参考基因
                    if ref_base_indel_next is None:
                        continue
                    ref_base_indel_next = ref_base_indel_next.lower()
                    pileup_list_indel_next = b.pileup_column(bam_file, chr, pos + 1, pos + 2)
                    if pileup_list_indel_next is None or pileup_list_indel_next[0] is None:
                        continue
                    pileup_list_indel_next = pileup_list_indel_next[0]
                    ## 如果跟参考基因组不同就用-1表示，如果相同就用0表示
                    indel = [-1 if i != ref_base_indel_next else 0 for i in pileup_list_indel_next]
                    indel_list = [indel[i] + indel_list[i] for i in range(min(len(indel), len(indel_list)))]
                    genotype_list = [ref_base_genotype + '-' if i < 0 else ref_base_genotype for i in indel_list]
                    genotype_list = [self.__str_to_int(i) for i in genotype_list]

                    ## 碱基改变系数
                    # base_coefficient = [1 if i < 0 else 0 for i in indel_list]
                    # base_coefficient = round(sum(base_coefficient) / len(indel_list), 3)


                ## padded_list, 对数据进行规整
                ref_var_list_padded = self.__padded_fill(ref_var_list, self.padded_maxlen)
                indel_list_padded = self.__padded_fill(indel_seq, self.padded_maxlen)
                genotype_list_padded = self.__padded_fill(genotype_list, self.padded_maxlen)


                var_ind_gen = ref_var_list_padded + indel_list_padded + genotype_list_padded
                # print(var_ind_gen)

                ## 处理label, 碱基label对应的和self.d里类似, indel label 用0，1，2, 基因型label用0，1，2，1，是基因型*/*的加和
                ### 处理碱基label (20种)
                indel_sum = sum(indel_list)
                if indel_sum != 0:  ## 如果是indel变异，则碱基标签为众数
                    ref_var_label = stats.mode(ref_var_list)[0][0] ## 找到众数作为标签, 众数不具有代表行
                elif indel_sum == 0: ## 如果是snp变异，碱基标签为vcf给出标签组合
                    ref_var_label = self.__str_to_int(str(ref_alts).lower() + str(alts).lower())

                ### 处理indel的label (3种)
                indel_sum = sum(indel_list)
                if indel_sum == 0:  ## 不是indel
                    indel_label = 0
                if indel_sum < 0:   ## 缺失
                    indel_label = 1
                if indel_sum > 0:   ## 插入
                    indel_label = 2

                ### 处理基因型label  (4种)
                if label == (0, 1):
                    g_label = 1  ## 杂合变异
                elif label == (1, 1):
                    g_label = 2  ## 纯合变异
                elif label == (1, 2):
                    g_label = 1  ## 杂合变异

                v_label = [ref_var_label, indel_label, g_label]
                # s_c_p = sample + '_' + chr + '_' + str(pos) + '_' + str(base_coefficient)
                s_c_p = sample + '_' + chr + '_' + str(pos)
                variant_sample = (s_c_p, var_ind_gen, v_label)
                samples_data.append(variant_sample)

            elif int(self.mode) == 2: ## 生成非变异数据
                ## 非变异数据
                pos = int(pos) + 1
                while pos:
                    normal_pileup_list = b.pileup_column(bam_file, chr, pos, pos + 1)
                    ref_base = fa.ref_atcg(fasta_file, chr, pos, pos + 1)
                    ## 如果不是正常序列就跳出
                    if normal_pileup_list is None or normal_pileup_list[0] is None or ref_base is None or normal_pileup_list[1] is None:
                        pos += 1
                        continue

                    ## 处理该位点的碱基序列
                    normal_seq = [item.lower() for item in normal_pileup_list[0]]
                    ref_base = ref_base.lower()
                    norepeat = set(normal_seq)  ## 去重
                    # print(norepeat)
                    if len(norepeat) == 1 and list(norepeat)[0] == ref_base:
                        ref_norm_list = [ref_base + i for i in normal_seq]
                        ref_norm_list = [self.__str_to_int(i) for i in ref_norm_list]

                        ## 处理该位点的indel序列
                        indel_norm_list = normal_pileup_list[1]
                        indel_norm_seq = [self.__indel_to_int(i) for i in indel_norm_list]
                        # print('normal', indel_norm_list)

                        ## 处理该位点的genotype序列
                        indel_norm_list_sum = sum(indel_norm_list)
                        if indel_norm_list_sum == 0:
                            # ref_norm_base_genotype = fa.ref_atcg(fasta_file, chr, pos, pos + 1)  ## 读下一个位置
                            # ref_norm_base_genotype = ref_norm_base_genotype.lower()
                            # genotype_norm_list = [ref_norm_base_genotype for i in indel_norm_list]
                            # genotype_norm_list = [self.__str_to_int(i, type=2) for i in genotype_norm_list]
                            genotype_norm_list = ref_norm_list

                            ## 碱基改变系数
                            # base_coefficient = round(indel_norm_list_sum / len(indel_norm_list), 3)

                        ref_norm_list_padded = self.__padded_fill(ref_norm_list, self.padded_maxlen)
                        indel_norm_list_padded = self.__padded_fill(indel_norm_seq, self.padded_maxlen)
                        genotype_norm_list_padded = self.__padded_fill(genotype_norm_list, self.padded_maxlen)

                        nor_ind_gen = ref_norm_list_padded + indel_norm_list_padded + genotype_norm_list_padded


                        ## 处理标签
                        ### 处理碱基label (20种), 如果没有变异，则众数就是组合
                        ref_norm_label = stats.mode(ref_norm_list)[0][0]  ## 找到众数作为标签
                        ### 处理indel
                        indel_norm_sum = sum(indel_norm_list)
                        if indel_norm_sum == 0:
                            indel_norm_label = 0
                        if indel_norm_sum < 0:
                            indel_norm_label = 1
                        if indel_norm_sum > 0:
                            indel_norm_label = 2
                        ### 处理基因型标签
                        g_norm_label = 0

                        n_label = [ref_norm_label, indel_norm_label, g_norm_label]

                        # normal_s_c_p = sample + '_' + chr + '_' + str(pos) + '_' + str(base_coefficient)
                        normal_s_c_p = sample + '_' + chr + '_' + str(pos)
                        normal_sample = (normal_s_c_p, nor_ind_gen, n_label)
                        samples_data.append(normal_sample)
                        break
                    pos += 1

        np.save(self.data_filename, samples_data)
        return samples_data


    ## 生成未知测试数据
    def test_pos(self):
        samples_data = []

        b = BAM()
        bam_file = b.readfile(self.bam_filename)

        fa = FASTA()
        fasta_file = fa.readfile(self.fasta_filename)

        t = TREGION()
        region_file = t.readfile(self.region_filename)
        region = t.region_info(region_file)

        for rec in region:
            # print('radom test')
            sample = rec[0]
            chr = rec[1]
            l_pos = rec[2]
            r_pos = rec[3]
            for pos in range(int(l_pos), int(r_pos) + 1):
                seq_list = b.pileup_column(bam_file, chr, pos, pos + 1)
                ref_base = fa.ref_atcg(fasta_file, chr, pos, pos + 1)

                if ref_base is None or seq_list is None or seq_list[0] is None or seq_list[1] is None:
                    continue

                ## 生成碱基序列
                seq_base_list = ['d' if item == '' else item for item in seq_list[0]]
                seq_lower_list = [item.lower() for item in seq_base_list]
                ref_base = ref_base.lower()
                ref_test_list = [ref_base + i for i in seq_lower_list]
                ref_test_list = [self.__str_to_int(i) for i in ref_test_list]
                ## 生成indel序列
                indel_test_list = seq_list[1]
                indel_test_seq = [self.__indel_to_int(i) for i in indel_test_list]

                ## 生成genotype序列
                indel_test_list_sum = sum(indel_test_list)
                ref_base_genotype_test = ref_base
                ref_base_genotype_test = ref_base_genotype_test.lower()
                if indel_test_list_sum == 0:
                    genotype_test_list = ref_test_list
                    ## 计算碱基改变系数
                    # base_coefficient = [1 if ref_base != i else 0 for i in seq_lower_list]
                    # base_coefficient = round(sum(base_coefficient) / len(base_coefficient), 3)

                    ## 计算indel跨越了多少个位点
                    indel_value = 0

                ## indel 基因型处理
                elif indel_test_list_sum > 0:
                    genotype_test_list = [ref_base_genotype_test + '+' if i > 0 else ref_base_genotype_test for i in indel_test_list]
                    genotype_test_list = [self.__str_to_int(i) for i in genotype_test_list]

                    ## 计算碱基改变系数
                    # base_coefficient = [1 if i > 0 else 0 for i in indel_test_list]
                    # base_coefficient = round(sum(base_coefficient) / len(indel_test_list), 3)

                    ## 计算indel跨越了多少个位点
                    indel_value = np.max(indel_test_list)

                elif indel_test_list_sum < 0:
                    ## 如果是缺失indel要特判一下
                    ref_base_indel_next_test = fa.ref_atcg(fasta_file, chr, pos + 1, pos + 2)  ## 取indel缺失位置位置
                    if ref_base_indel_next_test is None:
                        continue
                    ref_base_indel_next_test = ref_base_indel_next_test.lower()
                    pileup_list_indel_next_test = b.pileup_column(bam_file, chr, pos + 1, pos + 2)
                    if pileup_list_indel_next_test is None or pileup_list_indel_next_test[0] is None:
                        continue
                    pileup_list_indel_next_test = pileup_list_indel_next_test[0]
                    ## 如果跟参考基因组不同就用-1表示，如果相同就用0表示
                    indel = [-1 if i != ref_base_indel_next_test else 0 for i in pileup_list_indel_next_test]
                    indel_test_list = [indel[i] + indel_test_list[i] for i in range(min(len(indel), len(indel_test_list)))]

                    genotype_test_list = [ref_base_genotype_test + '-' if i < 0 else ref_base_genotype_test for i in indel_test_list]
                    genotype_test_list = [self.__str_to_int(i) for i in genotype_test_list]

                    ## 计算碱基改变系数
                    # base_coefficient = [1 if i < 0 else 0 for i in indel_test_list]
                    # base_coefficient = round(sum(base_coefficient) / len(indel_test_list), 3)

                    ## 计算indel跨越了多少个位点
                    indel_value = - np.min(indel_test_list)

                ref_test_list_padded = self.__padded_fill(ref_test_list, self.padded_maxlen)
                indel_test_list_padded = self.__padded_fill(indel_test_seq, self.padded_maxlen)
                genotyp_test_list_padded = self.__padded_fill(genotype_test_list, self.padded_maxlen)

                test_seq = ref_test_list_padded + indel_test_list_padded + genotyp_test_list_padded
                # s_c_p = sample + '_' + chr + '_' + str(pos) + '_' + str(base_coefficient) + '_' + str(seq_list[2]) + '_' + str(seq_list[3])
                s_c_p = sample + '_' + chr + '_' + str(pos) + '_' + str(indel_value) + '_' + str(seq_list[2]) + '_' + str(seq_list[3] )
                samples_data.append((s_c_p, test_seq))

        np.save(self.data_filename, samples_data)
        return samples_data

    def fetch_row(self, bam_file, chr_id, start, end, index, indel_value):
        i = 0
        for rec in bam_file.fetch(chr_id, start - 1 , end - 1):
            if i != index:
                i += 1
                continue

            seq = list(rec.seq)
            reference = rec.get_reference_sequence()
            reference = list(reference)
            pairs = rec.get_aligned_pairs()
            if indel_value > 0: ## 插入
                indel_insertion = ""
                for item in pairs:
                    if start in item and None not in item:
                        ref = reference[item[0]]   ##找到indel插入的参考基因
                        for i in range(indel_value):
                            indel_insertion += seq[item[0] + i + 1]  ## 找到后边插入的基因是什么
                        indel_insertion = ref + indel_insertion
                        re = ref.upper() + '-' + indel_insertion.upper()
                        return re

            elif indel_value < 0:
                indel_deletion = ""
                for item in pairs:
                    if start in item and None not in item:
                        ref = reference[item[0]]  ## 找到indel缺失的参考基因
                        for i in range(-indel_value):
                            indel_deletion += reference[item[0] + i + 1] ## 找到缺失的参考基因是什么

                        indel_deletion = ref + indel_deletion
                        re = indel_deletion.upper() + '-' + ref.upper()
                        return re


    def test_pos_all(self):
        samples_data = []

        b = BAM()
        bam_file = b.readfile(self.bam_filename)

        fa = FASTA()
        fasta_file = fa.readfile(self.fasta_filename)

        t = TREGION()
        region_file = t.readfile(self.region_filename)
        region = t.region_info(region_file)

        for rec in region:
            print('generation test data')
            sample = rec[0]
            chr = rec[1]
            l_pos = int(rec[2]) - 1
            r_pos = int(rec[3]) - 1
            for rec in bam_file.pileup(chr, l_pos, r_pos):
                print(rec.next().pos)
                try:
                    pos = rec.pos + 1  ## 参考基因位点
                    base_list = rec.get_query_sequences()
                    indel_list = [int(tmp.indel) for tmp in rec.pileups]

                    sum_indel_list = sum(indel_list)
                    if sum_indel_list == 0:  ## 如果是SNP
                        re = '0'
                    elif sum_indel_list < 0:
                        indel_index = np.argmin(indel_list)
                        indel_value = np.min(indel_list)
                        re = self.fetch_row(bam_file, chr, rec.pos, rec.pos + 1, indel_index, indel_value)
                    elif sum_indel_list > 0:
                        indel_value = np.max(indel_list)
                        indel_index = np.argmax(indel_list)
                        re = self.fetch_row(bam_file, chr, rec.pos, rec.pos + 1, indel_index, indel_value)

                    base_ad = Counter(indel_list)
                    ad = []  ## 计算不同等位基因数量
                    for k, v in base_ad.items():
                        if k != 0:
                            ad.append(v)
                    dp = len(indel_list)  ## 总的映射深度
                    d = dp - sum(ad)  ## 与参考基因相同的数量
                    if sum(indel_list) != 0:
                        ad = [str(i) for i in ad]
                        ad = ",".join(ad)  ## 每种等位基因的深度
                        ad = str(d) + ',' + ad
                        ad_dp = ad + '-' + str(dp)
                    else:
                        ad_dp = str(dp) + '-' + str(dp)

                    pileup_list = [base_list, indel_list, ad_dp, re]

                    ####################################################################
                    seq_list = pileup_list
                    ref_base = fa.ref_atcg(fasta_file, chr, pos, pos + 1)

                    if ref_base is None or seq_list is None or seq_list[0] is None or seq_list[1] is None:
                        continue

                    ## 生成碱基序列
                    seq_base_list = ['d' if item == '' else item for item in seq_list[0]]
                    seq_lower_list = [item.lower() for item in seq_base_list]
                    ref_base = ref_base.lower()
                    ref_test_list = [ref_base + i for i in seq_lower_list]
                    ref_test_list = [self.__str_to_int(i) for i in ref_test_list]
                    ## 生成indel序列
                    indel_test_list = seq_list[1]
                    indel_test_seq = [self.__indel_to_int(i) for i in indel_test_list]

                    ## 生成genotype序列
                    indel_test_list_sum = sum(indel_test_list)
                    ref_base_genotype_test = ref_base
                    ref_base_genotype_test = ref_base_genotype_test.lower()
                    if indel_test_list_sum == 0:
                        genotype_test_list = ref_test_list
                        ## 计算碱基改变系数
                        # base_coefficient = [1 if ref_base != i else 0 for i in seq_lower_list]
                        # base_coefficient = round(sum(base_coefficient) / len(base_coefficient), 3)

                        ## 计算indel跨越了多少个位点
                        indel_value = 0

                    ## indel 基因型处理
                    elif indel_test_list_sum > 0:
                        genotype_test_list = [ref_base_genotype_test + '+' if i > 0 else ref_base_genotype_test for i in
                                              indel_test_list]
                        genotype_test_list = [self.__str_to_int(i) for i in genotype_test_list]

                        ## 计算碱基改变系数
                        # base_coefficient = [1 if i > 0 else 0 for i in indel_test_list]
                        # base_coefficient = round(sum(base_coefficient) / len(indel_test_list), 3)

                        ## 计算indel跨越了多少个位点
                        indel_value = np.max(indel_test_list)

                    elif indel_test_list_sum < 0:
                        ## 如果是缺失indel要特判一下
                        ref_base_indel_next_test = fa.ref_atcg(fasta_file, chr, pos + 1, pos + 2)  ## 取indel缺失位置位置
                        if ref_base_indel_next_test is None:
                            continue
                        ref_base_indel_next_test = ref_base_indel_next_test.lower()
                        pileup_list_indel_next_test = b.pileup_column(bam_file, chr, pos + 1, pos + 2)
                        if pileup_list_indel_next_test is None or pileup_list_indel_next_test[0] is None:
                            continue
                        pileup_list_indel_next_test = pileup_list_indel_next_test[0]
                        ## 如果跟参考基因组不同就用-1表示，如果相同就用0表示
                        indel = [-1 if i != ref_base_indel_next_test else 0 for i in pileup_list_indel_next_test]
                        indel_test_list = [indel[i] + indel_test_list[i] for i in
                                           range(min(len(indel), len(indel_test_list)))]

                        genotype_test_list = [ref_base_genotype_test + '-' if i < 0 else ref_base_genotype_test for i in
                                              indel_test_list]
                        genotype_test_list = [self.__str_to_int(i) for i in genotype_test_list]

                        ## 计算碱基改变系数
                        # base_coefficient = [1 if i < 0 else 0 for i in indel_test_list]
                        # base_coefficient = round(sum(base_coefficient) / len(indel_test_list), 3)

                        ## 计算indel跨越了多少个位点
                        indel_value = - np.min(indel_test_list)

                    ref_test_list_padded = self.__padded_fill(ref_test_list, self.padded_maxlen)
                    indel_test_list_padded = self.__padded_fill(indel_test_seq, self.padded_maxlen)
                    genotyp_test_list_padded = self.__padded_fill(genotype_test_list, self.padded_maxlen)

                    test_seq = ref_test_list_padded + indel_test_list_padded + genotyp_test_list_padded
                    # s_c_p = sample + '_' + chr + '_' + str(pos) + '_' + str(base_coefficient) + '_' + str(seq_list[2]) + '_' + str(seq_list[3])
                    s_c_p = sample + '_' + chr + '_' + str(pos) + '_' + str(indel_value) + '_' + str(
                        seq_list[2]) + '_' + str(seq_list[3])

                    print((s_c_p, test_seq))
                    samples_data.append((s_c_p, test_seq))
                except Exception as e:
                    pass
                continue

            np.save(self.data_filename, samples_data)
            return samples_data

