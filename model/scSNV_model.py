from keras.models import Model
from keras.optimizers import Adam
from keras.layers import Input, Embedding, LSTM, Bidirectional, Dense, BatchNormalization
from keras.layers import Dropout
from keras.utils import multi_gpu_utils


class SCSNVMODEL:
    def __init__(self, input_shape=(78, ), n_labels=2, n_lstm_outdim=78, word_maxlen=78, em_inputdim=15, em_outdim=15, lstm_layers=1, dense_layers=1, dense_num=64, drop_out=0.2, lr=0.001, gpus=2):
        self.input_shape = input_shape
        self.n_labels = n_labels
        self.word_maxlen = word_maxlen
        self.n_lstm_outdim = n_lstm_outdim
        self.em_inputdim = em_inputdim
        self.em_outdim =em_outdim
        self.lstm_layers = lstm_layers
        self.dense_layers = dense_layers
        self.dense_num = dense_num
        self.drop_out = drop_out
        self.init_lr = lr
        self.gpus = gpus

    def model_construct(self):
        ## input_shape 代表句子的最大长度和input_length要对应
        inputs = Input(self.input_shape)
        '''
        input_dim: 词索引的种类，比如将所以词映射的索引是1-19之间，则input_dim = 20, 根据需要任意指定
        output_dim: embedding后每个词嵌入到一个多少维的向量中，比如将一个很高维的词索引嵌入到一个低维的索引中
        input_length: 句子最大长度，一个句子中最多能够有多少个词
        '''
        em = Embedding(input_dim=self.em_inputdim, output_dim=self.em_outdim, input_length=self.word_maxlen)(inputs)
        # print('em.shape', em.shape)
        bi_lstm = em
        for i in range(self.lstm_layers):
            if i == self.lstm_layers - 1:
                ## 最后一层只需要最后一个时刻的输入, 所以return_sequences=False
                bi_lstm = Bidirectional(LSTM(units=self.n_lstm_outdim*(i + 1), return_sequences=False))(bi_lstm)
                bi_lstm = BatchNormalization()(bi_lstm)
                # bi_lstm = LSTM(units=self.n_lstm_outdim*(i + 1), return_sequences=False)(bi_lstm)
                break

            ## return_sequences = True 不到最后一层, LSTM的下一层要用到上一层每个时刻的输入, 所以return_sequences设置为True
            bi_lstm = Bidirectional(LSTM(units=self.n_lstm_outdim*(i + 1), return_sequences=True))(bi_lstm)
            # bi_lstm = LSTM(units=self.n_lstm_outdim*(i + 1), return_sequences=True)(bi_lstm)
            bi_lstm = BatchNormalization()(bi_lstm)
            bi_lstm = Dropout(self.drop_out)(bi_lstm)

        dense = bi_lstm
        for i in range(self.dense_layers - 1, -1, -1):
            # print(i)
            dense = Dense(self.dense_num * (i + 1))(dense)

        outputs = Dense(self.n_labels, activation='softmax')(dense)
        model = Model(inputs, outputs)
        model.summary()
        # model = multi_gpu_utils(model, multi_gpu_utils=self.gpus)
        model.compile(optimizer=Adam(lr=self.init_lr), loss='categorical_crossentropy', metrics=['acc'])
        return model


if __name__ == '__main__':
    m = SCSNVMODEL()
    m.model_construct()






