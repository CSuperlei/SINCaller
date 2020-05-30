import numpy as np
import keras
from keras.preprocessing.text import one_hot
from keras.preprocessing.sequence import pad_sequences
from keras.utils.np_utils import to_categorical

class DataGenerator(keras.utils.Sequence):
    def __init__(self, samples_data, batch_size=64, shuffle=True, vocab_size=20, word_maxlen=78, label_len=4):
        self.samples_data = samples_data
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.sendin = []
        self.vocab_size = vocab_size
        self.word_maxlen = word_maxlen
        self.label_len = label_len
        self.indexes = None
        self.on_epoch_end()


    def __len__(self):
        ## 每一轮训练包含多少个batch
        return int(np.floor(len(self.samples_data) / self.batch_size))

    ## 每训练完打乱一次样本
    def on_epoch_end(self):
        ## 为每一个样本赋一个索引值
        self.indexes = np.arange(len(self.samples_data))
        if self.shuffle:
            ## 打乱所有样本的索引值
            np.random.shuffle(self.indexes)

    def __getitem__(self, index):
        # print('index', index)
        if (index + 1) * self.batch_size < len(self.samples_data):
            idx = self.indexes[index*self.batch_size : (index + 1)*self.batch_size]
            # print('idx', idx)
            ## 存放每个batch送入网络的索引值
            # self.sendin.append((idx))

            X, y = self.__data_generation(idx)

            return X, y

    def __data_generation(self, idx):
        ## 处理数据
        X = np.empty((self.batch_size, self.word_maxlen))
        y = np.empty((self.batch_size, self.label_len))
        for i, item in enumerate(idx):
            sample = self.samples_data[item]
            info = sample[0]
            self.sendin.append(info)
            ref = list(sample[1][0])
            seq = list(sample[1][1])
            data = ref + seq
            print('data', data)
            encoded_docs = [one_hot(d, self.vocab_size) for d in data]
            print('en', encoded_docs)
            padded_docs = pad_sequences(encoded_docs, maxlen=self.word_maxlen, padding='post')
            print('padded_docs shape', padded_docs.shape)
            label = sample[2]
            if label == (0, 0):
                label = 0
            elif label == (0, 1):
                label = 1
            elif label == (1, 1):
                label = 2
            elif label == (1, 2):
                label = 3

            label = to_categorical(label, num_classes=self.label_len)
            X[i, ] = padded_docs
            y[i, ] = label

        return X, y


    def get_sendin_content(self):
        return self.sendin