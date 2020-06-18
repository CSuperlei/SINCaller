import numpy as np
from keras.preprocessing.sequence import pad_sequences
from keras.utils.np_utils import to_categorical


class TEST:
    def __init__(self, samples_data, test_type=1, word_maxlen=234, label_base=20, label_indel=3, label_genotype=3):
        self.samples_data = samples_data
        self.word_maxlen = word_maxlen
        self.test_type = test_type
        self.label_base = label_base
        self.label_indel = label_indel
        self.label_genotype = label_genotype
        self.sendin = []

    def data_generator(self):
        if self.test_type == 1:  ## test_type == 1 生成有标签的数据，便于evaluate
            batch_data = []
            label_data1 = []
            label_data2 = []
            label_data3 = []
            for i, item in enumerate(self.samples_data):
                info = item[0]
                self.sendin.append(info)
                i_data = item[1]
                print(i_data)
                batch_data.append(i_data)
                label_base = item[2][0]
                label_indel = item[2][1]
                label_genotype = item[2][2]
                print(label_base)
                print(label_indel)
                print(label_genotype)
                label_data1.append(label_base)
                label_data2.append(label_indel)
                label_data3.append(label_genotype)

            padded_docs = pad_sequences(batch_data, maxlen=self.word_maxlen, padding='post')
            label_data1 = to_categorical(label_data1, num_classes=self.label_base)
            label_data2 = to_categorical(label_data2, num_classes=self.label_indel)
            label_data3 = to_categorical(label_data3, num_classes=self.label_genotype)
            X = np.array(padded_docs)
            y1 = np.array(label_data1)
            y2 = np.array(label_data2)
            y3 = np.array(label_data3)

            return X, {'outputs_base': y1, 'outputs_indel': y2, 'outputs_genotype': y3}

        elif self.test_type == 2:  ## 生成没有标签的数据
            batch_data = []
            for i, item in enumerate(self.samples_data):
                info = item[0]
                self.sendin.append(info)
                i_data = item[1]
                batch_data.append(i_data)

            padded_docs = pad_sequences(batch_data, maxlen=self.word_maxlen, padding='post')
            X = np.array(padded_docs)
            return X

    def get_sendin(self):
        return self.sendin




