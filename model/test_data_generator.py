import numpy as np
from keras.preprocessing.sequence import pad_sequences
from keras.utils.np_utils import to_categorical


class TEST:
    def __init__(self, samples_data, test_type=1, word_maxlen=78, label_len=2):
        self.samples_data = samples_data
        self.word_maxlen = word_maxlen
        self.test_type = test_type
        self.label_len = label_len
        self.d = {'aa': 1, 'at': 2, 'ac': 3, 'ag': 4, 'ad': -1,
                  'tt': 5, 'ta': 6, 'tc': 7, 'tg': 8, 'td': -1,
                  'cc': 9, 'ca': 10, 'ct': 11, 'cg': 12, 'cd': -1,
                  'gg': 13, 'ga': 14, 'gc': 15, 'gt': 16, 'gd': -1,
                  }
        self.sendin = []

    def __str_to_int(self, s):
        r = self.d[s]
        return r

    def data_generator(self):
        if self.test_type == 1:
            test_data = []
            test_label = []
            for i, item in enumerate(self.samples_data):
                info = item[0]
                self.sendin.append(info)
                ref = item[1][0]
                seq = item[1][1]
                label = item[2]
                data = [ref + i for i in seq]
                i_data = [self.__str_to_int(i) for i in data]
                test_data.append(i_data)
                if label == (0, 0):
                    label = 0
                    test_label.append(label)
                elif label == (0, 1):
                    label = 1
                    test_label.append(label)
                elif label == (1, 1):
                    label = 1
                    test_label.append(label)
                elif label == (1, 2):
                    label = 1
                    test_label.append(label)

            padded_docs = pad_sequences(test_data, maxlen=self.word_maxlen, padding='post')
            label_data = to_categorical(test_label, num_classes=self.label_len)
            X = np.array(padded_docs)
            y = np.array(label_data)
            return X, y

        elif self.test_type == 2:
            test_data = []
            for i, item in enumerate(self.samples_data):
                info = item[0]
                self.sendin.append(info)
                ref = item[1][0]
                seq = item[1][1]
                data = [ref + i for i in seq]
                i_data = [self.__str_to_int(i) for i in data]
                test_data.append(i_data)

            padded_docs = pad_sequences(test_data, maxlen=self.word_maxlen, padding='post')
            X = np.array(padded_docs)
            return X

    def get_sendin(self):
        return self.sendin




