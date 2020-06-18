from .train_data_generator import DataGenerator
from .scSNV_model import SCSNVMODEL
import time
from tensorflow.keras.callbacks import TensorBoard, EarlyStopping, ModelCheckpoint, LearningRateScheduler
import tensorflow.keras.backend as K


def training(samples_train_data, samples_val_data, epochs=30, generator_params=None, model_params=None, hdf5_file = False, hdf5_fliename=None, mcheckpoint_dir='/home/cailei/bio_project/nbCNV/train_log/model_checkpoint/', mtensorboard_dir='/home/cailei/bio_project/nbCNV/train_log/tensorboard_logs/'):
    print('sample_train_data', len(samples_train_data))
    print('sample_val_data', len(samples_val_data))
    training_generator = DataGenerator(samples_train_data)
    validation_generator = DataGenerator(samples_val_data)
    if generator_params:
        training_generator = DataGenerator(samples_train_data, **generator_params)
        validation_generator = DataGenerator(samples_val_data, **generator_params)

    m = SCSNVMODEL()
    if model_params:
        m = SCSNVMODEL(**model_params)

    model = m.model_construct()
    if hdf5_file and hdf5_fliename != None:
        model.load_weights(hdf5_fliename)

    cb_1 = EarlyStopping(min_delta=0, patience=30, verbose=0, mode='auto')
    cb_2 = ModelCheckpoint(filepath=mcheckpoint_dir+'model_{epoch:02d}.hdf5', verbose=0,
                                           save_best_only=False, save_weights_only=False, mode='auto', period=1)
    model_name = 'tensorboard_scSNV_{}'.format(int(time.time()))
    cb_3 = TensorBoard(log_dir=mtensorboard_dir+'{}'.format(model_name))

    def scheduler(epoch):
        # 每隔100个epoch，学习率减小为原来的1/10
        if epoch % 10 == 0 and epoch != 0:
            lr = K.get_value(model.optimizer.lr)
            print('current lr', lr)
            K.set_value(model.optimizer.lr, lr * 0.1)
            print("lr changed to {}".format(lr * 0.1))
        return K.get_value(model.optimizer.lr)

    cb_4 = LearningRateScheduler(scheduler)
    results = model.fit_generator(generator=training_generator,
                        validation_data=validation_generator,
                        epochs=epochs,
                        callbacks=[cb_1, cb_2, cb_3, cb_4])
    return training_generator.get_sendin_content()


