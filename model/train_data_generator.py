import numpy as np
import keras
from keras.preprocessing.text import one_hot
from keras.preprocessing.sequence import pad_sequences
from keras.utils.np_utils import to_categorical

class DataGenerator(keras.utils.Sequence):
    def __init__(self, samples_data, batch_size=128, shuffle=True, word_maxlen=78, label_base=2, label_indel=3, label_genotype=3):
        self.samples_data = samples_data
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.sendin = []
        self.word_maxlen = word_maxlen
        self.label_base = label_base
        self.label_indel = label_indel
        self.label_genotype = label_genotype
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
        # X = np.empty((self.batch_size, self.word_maxlen))
        # y = np.empty((self.batch_size, self.label_len))
        batch_data = []
        label_data1 = []
        label_data2 = []
        label_data3 = []
        for i, item in enumerate(idx):
            sample = self.samples_data[item]
            info = sample[0]
            self.sendin.append(info)
            i_data = sample[1]
            batch_data.append(i_data)
            label_base = sample[2]
            # label_base = sample[2][0]
            # label_indel = sample[2][1]
            # label_genotype = sample[2][2]
            label_data1.append(label_base)
            # label_data2.append(label_indel)
            # label_data3.append(label_genotype)

        # print(batch_data[-1])
        padded_docs = pad_sequences(batch_data, maxlen=self.word_maxlen, padding='post')
        # print('padded_docs shape', padded_docs.shape)
        print('padded_docs', padded_docs)
        label_data1 = to_categorical(label_data1, num_classes=self.label_base)
        # label_data2 = to_categorical(label_data2, num_classes=self.label_indel)
        # label_data3 = to_categorical(label_data3, num_classes=self.label_genotype)
        print('label_data1', label_data1)
        # print('label_data2', label_data2)
        # print('label_data3', label_data3)
        X = np.array(padded_docs)
        y1 = np.array(label_data1)
        # print(y1[-1])
        # y2 = np.array(label_data2)
        # y3 = np.array(label_data3)
        # return X, {'outputs_base': y1, 'outputs_indel': y2, 'outputs_genotype':y3}
        return X, {'outputs_base': y1}

    def get_sendin_content(self):
        return self.sendin