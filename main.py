import os
import argparse
import numpy as np
from data_process.data_preprocess import DATAPROCESS
from model.train_process import training
from model.test_process import testing
from keras.preprocessing.sequence import pad_sequences
from keras.utils.np_utils import to_categorical


def args_func():
    parser = argparse.ArgumentParser(description="scSNV software")
    parser.add_argument('--vcf', '-v', help='vcf filename')
    parser.add_argument('--bam', '-b', help='bam filename')
    parser.add_argument('--fasta', '-fa', help='fasta, filename')
    parser.add_argument('--fastq', '-fq', help='fastq filename')
    parser.add_argument('--gpus', '-g', help='gpu number')
    parser.add_argument('--log', '-lo', help='log level')
    parser.add_argument('--data', '-d', help='data filename')   ## 生成数据的名字
    parser.add_argument('--load', '-ld', help='load filename')  ## 加载数据
    parser.add_argument('--test', '-tm', help='test mode 1 is generator_test; mode 2 is batch test ')
    parser.add_argument('--mode', '-m', help='mode 1 is training; mode 2 is tesing; mode 3 is generate data', required=True)
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
    # b.pileup_column(bam_file, 'chr1', 2956920, 2956921)

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
            sendin = training(samples_train_data, samples_val_data)
            # print(sendin)  ## 打印训练数据
        else:
            print('data is empty')
            return

    elif int(mode) == 2:  ## 测试
        load_filename = args.load
        if load_filename is not None:
            samples_data = np.load(load_filename, allow_pickle=True)
       
        if samples_data is not None and int(args.test) == 1:
            generator_params = { 'shuffle': False }
            sendin = testing(samples_data, test_model=1, generator_params=generator_params)

        elif int(args.test) == 2:
            d = {'aa': 1, 'at': 2, 'ac': 3, 'ag': 4,
                 'tt': 5, 'ta': 6, 'tc': 7, 'tg': 8,
                 'cc': 9, 'ca': 10, 'ct': 11, 'cg': 12,
                 'gg': 13, 'ga': 14, 'gc': 15, 'gt': 16}
            docs = ['at at',
                    'aa aa',
                    'gt gt gt gt gt',
                    'cc cc cc cc cc',
                    'ag ag',
                    'tt tt',
                    'tc tc',
                    'cc cc',
                    'gc gc gc',
                    'aa aa aa',
                    'gg gg gt',
                    'tt tt tt',
                    'ag ag ag ag ag ag ag ag ag ag ag ag ag ag',
                    'gg gg gg gg gg gg gg gg gg gg gg gg gg gg',
                    'ga ga ga',
                    'cc cc cc',
                    ]

            for i, item in enumerate(docs):
                t = []
                for tmp in item.split(' '):
                    r = d[tmp]
                    t.append(r)
                docs[i] = t
                t = []

            labels = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
            # labels = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
            ohl = to_categorical(labels, num_classes=2)
            print('ohl', ohl)
            max_length = 20
            padded_docs = pad_sequences(docs, maxlen=max_length, padding='post')
            print('padded_docs', padded_docs)
            print('padded_docs shape', padded_docs.shape)
            sample_test_data = [padded_docs, ohl]
            # print(sample_test_data)
            testing(sample_test_data, test_model=2)

    elif int(mode) == 3:  ## 生成数据
        vcf_filename = args.vcf
        bam_filename = args.bam
        fasta_filename = args.fasta
        data_filename = args.data
        if data_filename is not None:
            d = DATAPROCESS(vcf_filename, bam_filename, fasta_filename, data_filename)
            samples_data = d.dataproc()
            print(samples_data)


if __name__ == '__main__':
    main()



