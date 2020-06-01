import sys
# import pysam
# from pysam import VariantFile
import argparse
from keras.utils.np_utils import to_categorical
from keras.preprocessing.text import one_hot
from keras.preprocessing.sequence import pad_sequences

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

    print(docs)
    labels = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
    # labels = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    ohl = to_categorical(labels, num_classes=2)
    print('ohl', ohl)
    # integer encode the documents
    vocab_size = 17
    # encoded_docs = [one_hot(d, vocab_size) for d in docs]
    # print('en', encoded_docs)
    # pad documents to a max length of 4 words
    max_length = 20
    padded_docs = pad_sequences(docs, maxlen=max_length, padding='post')
    print('padded_docs', padded_docs)
    print('padded_docs shape', padded_docs.shape)


if __name__ == '__main__':
    main()