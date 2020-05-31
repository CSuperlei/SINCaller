import keras
import time
from .data_generator import DataGenerator
from .scSNV_model import SCSNVMODEL


def training(samples_train_data, samples_val_data, epochs=20, generator_params=None, model_params=None, hdf5_file = False, hdf5_fliename=None, mcheckpoint_dir='/home/cailei/bio_project/nbCNV/train_log/model_checkpoint/', mtensorboard_dir='./tensorboard_logs/'):
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

    cb_1 = keras.callbacks.EarlyStopping(min_delta=0, patience=30, verbose=0, mode='auto')
    cb_2 = keras.callbacks.ModelCheckpoint(filepath=mcheckpoint_dir+'model_{epoch:02d}.hdf5', verbose=0,
                                           save_best_only=False, save_weights_only=False, mode='auto', period=1)
    # model_name = 'tensorboard_scSNV_{}'.format(int(time.time()))
    # cb_3 = keras.callbacks.TensorBoard(log_dir=mtensorboard_dir+'{}'.format(model_name))
    results = model.fit_generator(generator=training_generator,
                                  validation_data=validation_generator,
                                  epochs=epochs,
                                  nb_worker=1,
                                  callbacks=[cb_1, cb_2 ])

    return training_generator.get_sendin_content()


