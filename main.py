import os
import argparse
import numpy as np
from data_process.data_preprocess import DATAPROCESS
from model.train_process import training
from model.test_process import testing
from keras.preprocessing.text import one_hot
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
    parser.add_argument('--data', '-d', help='data filename')
    parser.add_argument('--load', '-ld', help='load filename')
    parser.add_argument('--mode', '-m', help='mode 1 is training; mode 2 is tesing', required=True)
    args = parser.parse_args()
    return args


def main():
    args = args_func()

    # vcf_filename = args.vcf
    # v = VCF()
    # v.readfile(vcf_filename)

    # bam_filename = args.bam
    #     # b = BAM()
    #     # bam_file = b.readfile(bam_filename)
    #     # b.pileup_column(bam_file, 'chr1', 2956920, 2956921)

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

    if int(mode) == 1:
        vcf_filename = args.vcf
        bam_filename = args.bam
        fasta_filename = args.fasta
        data_filename = args.data
        load_filename = args.load

        if load_filename is not None:
            samples_data = np.load(load_filename, allow_pickle=True)
        elif data_filename is not None:
            d = DATAPROCESS(vcf_filename, bam_filename, fasta_filename, data_filename)
            samples_data = d.dataproc()

        # print('samples_data shape', samples_data.shape)

        if samples_data is not None:
            # print(samples_data)
            sep = int(np.floor(len(samples_data) * 0.8))
            samples_train_data = samples_data[:sep]
            samples_val_data = samples_data[sep:]
            sendin = training(samples_train_data, samples_val_data)
            # print(sendin)
        else:
            print('data is empty')

    elif int(mode) == 2:

        docs = ['a, t, t',
                'a, a, a',
                'g, t, t, t, t, t',
                'c, c, c, c, c, c',
                'a, g, g',
                't, t, t',
                't, c, c',
                'c, c, c',
                'g, c, c, c',
                'a, a, a, a',
                'g, g, g, t',
                't, t, t, t',
                'a, g, g, g, g, g, g, g, g, g, g, g, g, g, g',
                'g, g, g, g, g, g, g, g, g, g, g, g, g, g, g',
                'g a a a',
                'c c c c',
                ]
        labels = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
        # labels = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        ohl = to_categorical(labels, num_classes=2)
        # print(ohl)
        # integer encode the documents
        vocab_size = 20
        encoded_docs = [one_hot(d, vocab_size) for d in docs]
        # print('en', encoded_docs)
        # pad documents to a max length of 4 words
        max_length = 78
        padded_docs = pad_sequences(encoded_docs, maxlen=max_length, padding='post')
        # print(padded_docs)
        # print('padded_docs shape', padded_docs.shape)
        sample_test_data = [padded_docs, ohl]
        print(sample_test_data)
        testing(sample_test_data)


if __name__ == '__main__':
    main()



