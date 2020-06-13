import os
import argparse
import numpy as np
from data_process.data_preprocess import DATAPROCESS
from data_process.data_combine import DATACOMBINE
from model.train_process import training
from model.test_process import testing

from keras.preprocessing.sequence import pad_sequences
from keras.utils.np_utils import to_categorical
from bam.bam_process import BAM


def args_func():
    parser = argparse.ArgumentParser(description="scSNV software")
    parser.add_argument('--vcf', '-v', help='vcf filename')
    parser.add_argument('--bam', '-b', help='bam filename')
    parser.add_argument('--fasta', '-fa', help='fasta, filename')
    parser.add_argument('--fastq', '-fq', help='fastq filename')
    parser.add_argument('--gpus', '-g', help='gpu number')
    parser.add_argument('--log', '-lo', help='log level')
    parser.add_argument('--data', '-d', help='data filename')   ## 生成数据的名字
    parser.add_argument('--data_model', '-dm', help='data generator mode; mode 1 gernerates variant data, mode 2 generates normal data')
    parser.add_argument('--load', '-ld', help='load filename')  ## 加载数据
    parser.add_argument('--region', '-r', help='region test filename')  ## 加载测试区域数据
    parser.add_argument('--test', '-tm', help='test mode 1 is generator_test; mode 2 is batch test; mode 3 is random data')
    parser.add_argument('--dc_origin', '-dco', help='data combine orgin') ## 合并不同标签数据
    parser.add_argument('--dc_target', '-dct', help='data combine target') ## 生成不同标签数据
    parser.add_argument('--mode', '-m', help='mode 1 is training; mode 2 is tesing; mode 3 is generate data; mode 4 is combine data', required=True)
    args = parser.parse_args()
    return args


def main():
    args = args_func()
    # vcf_filename = args.vcf
    # v = VCF()
    # v.readfile(vcf_filename)

    # bam_filename = args.bam
    # b = BAM()
    # bam_file = b.readfile(bam_filename)
    # b.pileup_column(bam_file, 'chr1', 54031886, 54031887)

    # fasta_filename = args.fasta
    # a = FASTA()
    # fasta_file = a.readfile(fasta_filename)
    # a.ref_atcg(fasta_file, 'chr1', 2956920, 2957922)

    # fastq_filename = args.fastq
    # q = FASTQ()
    # fastq_file = q.readfile(fastq_filename)
    # q.seq_atcg(fastq_file)

    gpus = args.gpus
    if gpus is not None:
        gpus = ",".join(list(str(gpus)))
        # print('gpus', gpus)
        os.environ["CUDA_VISIBLE_DEVICES"] = gpus

    log_level = args.log
    if log_level is not None:
        # print('log_level', log_level)
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = log_level

    mode = args.mode
    if int(mode) == 1:   ## 训练
        load_filename = args.load
        if load_filename is not None:
            samples_data = np.load(load_filename, allow_pickle=True)
        else:
            print('load filename is empty')
            return
        if samples_data is not None:
            sep = int(np.floor(len(samples_data) * 0.8))
            samples_train_data = samples_data[:sep]
            samples_val_data = samples_data[sep:]
            model_params = {
                'alpha_base': [[1],[3],[3],[3],[5],[1],[3],[3],[3],[5],[1],[3],[3],[3],[5],[1],[3],[3],[3],[5]],
                'alpha_indel': [[1],[5],[5]],
                'alpha_genotype': [[1],[5],[3],[7]]
            }
            sendin = training(samples_train_data, samples_val_data, model_params=model_params)
            # print(sendin)  ## 打印训练数据
        else:
            print('data is empty')
            return

    elif int(mode) == 2:  ## 测试
        if int(args.test) == 1:
            load_filename = args.load
            if load_filename is not None:
                samples_data = np.load(load_filename, allow_pickle=True)
            else:
                print('load filename is empyt')
                return
            if samples_data is not None and int(args.test) == 1:
                generator_params = { 'shuffle': False }
                sendin = testing(samples_data, test_model=1, generator_params=generator_params)

        elif int(args.test) == 2:
            # d = {
            #     'aa': 1, 'at': 2, 'ac': 3, 'ag': 4, 'ad': -1,
            #     'tt': 5, 'ta': 6, 'tc': 7, 'tg': 8, 'td': -1,
            #     'cc': 9, 'ca': 10, 'ct': 11, 'cg': 12, 'cd': -1,
            #     'gg': 13, 'ga': 14, 'gc': 15, 'gt': 16, 'gd': -1,
            # }
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

            samples_data = [
                ('SRR052047_chr2_49671216', ('a', ('t', 't')), (1, 1)),
                ('SRR052047_chr2_49671217', ('a', ('a', 'a')), (0, 0)),
                ('SRR052047_chr2_49671218', ('g', ('t', 't', 't', 't', 't')), (1, 1)),
                ('SRR052047_chr2_49671219', ('c', ('c', 'c', 'c', 'c', 'c')), (0, 0)),
                ('SRR052047_chr2_49671220', ('a', ('g', 'g')), (1, 1)),
                ('SRR052047_chr2_49671221', ('t', ('t', 't')), (0, 0)),
                ('SRR052047_chr2_49671222', ('t', ('c', 'c')), (1, 1)),
                ('SRR052047_chr2_49671223', ('c', ('c', 'c')), (0, 0)),
                ('SRR052047_chr2_49671224', ('g', ('c', 'c')), (1, 1)),
                ('SRR052047_chr2_49671225', ('a', ('a', 'a')), (0, 0)),
                ('SRR052047_chr2_49671226', ('g', ('g', 't')), (1, 1)),
                ('SRR052047_chr2_49671227', ('t', ('t', 't')), (0, 0)),
                ('SRR052047_chr2_49671228', ('a', ('g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g')), (1, 1)),
                ('SRR052047_chr2_49671229', ('g', ('g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g')), (0, 0)),
                ('SRR052047_chr2_49671230', ('g', ('a', 'a')), (1, 1)),
                ('SRR052047_chr2_49671231', ('c', ('c', 'c')), (0, 0)),
            ]
            # for i, item in enumerate(docs):
            #     t = []
            #     for tmp in item.split(' '):
            #         r = d[tmp]
            #         t.append(r)
            #     docs[i] = t
            #     t = []
            #
            # labels = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
            # # labels = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
            # ohl = to_categorical(labels, num_classes=2)
            # print('ohl', ohl)
            # max_length = 78
            # padded_docs = pad_sequences(docs, maxlen=max_length, padding='post')
            # print('padded_docs', padded_docs)
            # print('padded_docs shape', padded_docs.shape)
            # sample_test_data = [padded_docs, ohl]
            # print(sample_test_data)
            testing(samples_data, test_model=2)

        elif int(args.test) == 3:
            load_filename = args.load
            if load_filename is not None:
                samples_data = np.load(load_filename, allow_pickle=True)
            else:
                print('load filename is empyt')
                return
            testing(samples_data, test_model=3)


    elif int(mode) == 3:  ## 生成数据
        vcf_filename = args.vcf
        bam_filename = args.bam
        fasta_filename = args.fasta
        data_filename = args.data
        region_filename = args.region
        data_model = args.data_model
        ## 生成已知label数据
        if data_filename is not None and region_filename is None:
            d = DATAPROCESS(data_model, vcf_filename, bam_filename, fasta_filename, data_filename)
            samples_data = d.dataproc()
            print(samples_data)
        ## 生成未知label数据
        elif data_filename is not None and region_filename is not None:
            d = DATAPROCESS(data_model, vcf_filename=None, bam_filename=bam_filename, fasta_filename=fasta_filename, data_filename=data_filename, region_filename=region_filename)
            samples_data = d.test_pos()
            print(samples_data)

    elif int(mode) == 4:
        origin_dir = args.dc_origin
        targe_filename = args.dc_target
        dc = DATACOMBINE(origin_dir, targe_filename)
        dc.data_combine()


if __name__ == '__main__':
    main()



