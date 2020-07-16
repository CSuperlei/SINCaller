import tensorflow as tf
from tensorflow.keras.layers import *
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam, RMSprop
from tensorflow.keras.utils import plot_model, multi_gpu_model
from tensorflow.keras.layers import Embedding

class SCSNVMODEL:
    def __init__(self, input_shape=(3000, ), n_base_labels=20, n_indel_labels=3, n_genotype_labels=3, n_lstm_outdim=64, word_maxlen=3000, em_inputdim=36, em_outdim=36, lstm_layers=3, dense_layers=2, dense_num=64, drop_out=0.5, lr=0.001, gpus=2, alpha_base=[[1],[3],[3],[3],[5],[1],[3],[3],[3],[5],[1],[3],[3],[3],[5],[1],[3],[3],[3],[5]], alpha_indel=[[1],[5],[5]], alpha_genotype=[[1],[5],[3],[7]], gamma=2.0):
        self.input_shape = input_shape
        self.n_base_labels = n_base_labels
        self.n_indel_labels = n_indel_labels
        self.n_genotype_labels = n_genotype_labels
        self.word_maxlen = word_maxlen
        self.n_lstm_outdim = n_lstm_outdim
        self.em_inputdim = em_inputdim  ## 和one-hot的vocab_size可以相同
        self.em_outdim =em_outdim
        self.lstm_layers = lstm_layers
        self.dense_layers = dense_layers
        self.dense_num = dense_num
        self.drop_out = drop_out
        self.init_lr = lr
        self.gpus = gpus
        self.alpha_base = alpha_base
        self.alpha_indel = alpha_indel
        self.alpha_genotype = alpha_genotype
        self.gamma = gamma

    def multi_category_focal_loss1(self, alpha, gamma=2.0):
        """
        focal loss for multi category of multi label problem
        适用于多分类或多标签问题的focal loss
        alpha用于指定不同类别/标签的权重，数组大小需要与类别个数一致
        当你的数据集不同类别/标签之间存在偏斜，可以尝试适用本函数作为loss
        Usage:
         model.compile(loss=[multi_category_focal_loss1(alpha=[1,2,3,2], gamma=2)], metrics=["accuracy"], optimizer=adam)
        """
        epsilon = 1.e-7
        alpha = tf.constant(alpha, dtype=tf.float32)
        # alpha = tf.constant([[1],[1],[1],[1],[1]], dtype=tf.float32)
        # alpha = tf.constant_initializer(alpha)
        gamma = float(gamma)

        def multi_category_focal_loss1_fixed(y_true, y_pred):
            y_true = tf.cast(y_true, tf.float32)
            y_pred = tf.clip_by_value(y_pred, epsilon, 1. - epsilon)
            y_t = tf.multiply(y_true, y_pred) + tf.multiply(1 - y_true, 1 - y_pred)
            # ce = -tf.math.log(y_t)
            ce = -tf.log(y_t)
            weight = tf.pow(tf.subtract(1., y_t), gamma)
            fl = tf.matmul(tf.multiply(weight, ce), alpha)
            loss = tf.reduce_mean(fl)
            return loss

        return multi_category_focal_loss1_fixed

    def model_construct(self):
        ## input_shape 代表句子的最大长度和input_length要对应
        inputs = Input(self.input_shape)
        '''
        input_dim: 词索引的种类，比如将所以词映射的索引是1-19之间，则input_dim = 20, 根据需要任意指定
        output_dim: embedding后每个词嵌入到一个多少维的向量中，比如将一个很高维的词索引嵌入到一个低维的索引中
        input_length: 句子最大长度，一个句子中最多能够有多少个词
        '''
        em = Embedding(input_dim=self.em_inputdim, output_dim=self.em_outdim, input_length=self.word_maxlen)(inputs)

        split = Lambda(tf.split, arguments={'axis': 1, 'num_or_size_splits': 3})(em)
        lstm_base = split[0]
        lstm_indel = split[1]
        lstm_genotype = split[2]

        for i in range(self.lstm_layers):
            if i == self.lstm_layers - 1:
                ## 最后一层只需要最后一个时刻的输入, 所以return_sequences=False
                lstm_base = Bidirectional(LSTM(units=self.n_lstm_outdim*(i + 1), return_sequences=False))(lstm_base)
                lstm_base = BatchNormalization()(lstm_base)
                break

            ## return_sequences = True 不到最后一层, LSTM的下一层要用到上一层每个时刻的输入, 所以return_sequences设置为True
            lstm_base = Bidirectional(LSTM(units=self.n_lstm_outdim*(i + 1), return_sequences=True))(lstm_base)
            lstm_base = BatchNormalization()(lstm_base)
            lstm_base = Dropout(self.drop_out)(lstm_base)

        for i in range(self.lstm_layers):
            if i == self.lstm_layers - 1:
                ## 最后一层只需要最后一个时刻的输入, 所以return_sequences=False
                lstm_indel = Bidirectional(LSTM(units=self.n_lstm_outdim * (i + 1), return_sequences=False))(lstm_indel)
                lstm_indel = BatchNormalization()(lstm_indel)
                break

            ## return_sequences = True 不到最后一层, LSTM的下一层要用到上一层每个时刻的输入, 所以return_sequences设置为True
            lstm_indel = Bidirectional(LSTM(units=self.n_lstm_outdim * (i + 1), return_sequences=True))(lstm_indel)
            lstm_indel = BatchNormalization()(lstm_indel)
            lstm_indel = Dropout(self.drop_out)(lstm_indel)

        for i in range(self.lstm_layers):
            if i == self.lstm_layers - 1:
                ## 最后一层只需要最后一个时刻的输入, 所以return_sequences=False
                lstm_genotype = Bidirectional(LSTM(units=self.n_lstm_outdim * (i + 1), return_sequences=False))(lstm_genotype)
                lstm_genotype = BatchNormalization()(lstm_genotype)
                break

            ## return_sequences = True 不到最后一层, LSTM的下一层要用到上一层每个时刻的输入, 所以return_sequences设置为True
            lstm_genotype = Bidirectional(LSTM(units=self.n_lstm_outdim * (i + 1), return_sequences=True))(lstm_genotype)
            lstm_genotype = BatchNormalization()(lstm_genotype)
            lstm_genotype = Dropout(self.drop_out)(lstm_genotype)


        dense_base = lstm_base
        dense_indel = lstm_indel
        dense_genotype = lstm_genotype

        for i in range(self.dense_layers - 1, -1, -1):
            dense_base = Dense(self.dense_num * (i + 1))(dense_base)

        outputs_base = Dense(self.n_base_labels, activation='softmax', name='outputs_base')(dense_base)

        for i in range(self.dense_layers - 1, -1, -1):
            dense_indel = Dense(self.dense_num * (i + 1))(dense_indel)

        outputs_indel = Dense(self.n_indel_labels, activation='softmax', name='outputs_indel')(dense_indel)

        for i in range(self.dense_layers - 1, -1, -1):
            dense_genotype = Dense(self.dense_num * (i + 1))(dense_genotype)

        outputs_genotype = Dense(self.n_genotype_labels, activation='softmax', name='outputs_genotype')(dense_genotype)

        model = Model(inputs=inputs, outputs=[outputs_base, outputs_indel, outputs_genotype], name='model')
        ## 必须放到model.comile前边
        plot_model(model, expand_nested=True, show_shapes=True, show_layer_names=True, to_file='model_strcut.png', dpi=300)

        # model = multi_gpu_utils(model, multi_gpu_utils=self.gpus)
        # model.compile(optimizer=RMSprop(lr=self.init_lr), loss={'outputs_base': self.multi_category_focal_loss1(self.alpha_base, self.gamma), 'outputs_indel': self.multi_category_focal_loss1(self.alpha_indel, self.gamma), 'outputs_genotype':self.multi_category_focal_loss1(self.alpha_genotype, self.gamma)}, metrics=['acc'])
        model.compile(optimizer=Adam(lr=self.init_lr), loss={'outputs_base':'categorical_crossentropy', 'outputs_indel': 'categorical_crossentropy', 'outputs_genotype':'categorical_crossentropy'}, metrics=['acc'])
        model.summary()
        return model


if __name__ == '__main__':
    m = SCSNVMODEL()
    m.model_construct()


