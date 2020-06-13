import os
from keras import backend as K
import tensorflow as tf
import numpy as np
from keras.utils.np_utils import to_categorical

os.environ["CUDA_VISIBLE_DEVICES"] = '0'

def multi_category_focal_loss1(y_true, y_pred):
    epsilon = 1.e-7
    gamma = 2.0
    # alpha = tf.constant([[1],[1],[1],[1],[1]], dtype=tf.float32)
    alpha = tf.constant([[2],[1],[1],[1],[1],[2],[1],[1],[1],[1]], dtype=tf.float32)
    # alpha = tf.constant([[1],[1],[1],[1],[1]], dtype=tf.float32)

    y_true = tf.cast(y_true, tf.float32)
    y_pred = tf.clip_by_value(y_pred, epsilon, 1. - epsilon)
    y_t = tf.multiply(y_true, y_pred) + tf.multiply(1-y_true, 1-y_pred)
    ce = -tf.log(y_t)
    weight = tf.pow(tf.subtract(1., y_t), gamma)
    fl = tf.matmul(tf.multiply(weight, ce), alpha)
    loss = tf.reduce_mean(fl)
    return loss

def one_hot_test():
    a = [31, 32, 33, 34, 35,
         36, 37, 38, 39, 40]
    a = np.array(a)
    a = a - 31
    print(a)
    b = to_categorical(a, num_classes=10)
    print(b)

def numpy_con():
    a = [(1,2),(3,4),(5,6)]
    b = [(7,8),(8,9)]
    c = a + b
    print(c)
    a = np.array(a)
    b = np.array(b)
    c = np.concatenate((a, b), axis=0)
    c = list(c)
    print(c)

def main():

    # Y_true = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    # Y_pred = np.array([[0.93, 0.99, 0.8, 0.97, 0.15, 0.93, 0.99, 0.8, 0.97, 0.15], [0.1, 0.05, 0.1, 0.09, 0.9, 0.1, 0.05, 0.1, 0.09, 0.9]], dtype=np.float32)
    # print(K.eval(multi_category_focal_loss1(Y_true, Y_pred)))
    # one_hot_test()
    numpy_con()

if __name__ == '__main__':
    main()