from .scSNV_model import SCSNVMODEL
from .data_generator import DataGenerator

def testing(samples_test_data, test_model=1, model_params=None,  generator_params=None, hdf5_file=True, hdf5_fliename='/home/cailei/bio_project/nbCNV/train_log/model_checkpoint/model_12.hdf5'):
    print('samples test data', len(samples_test_data))
    m = SCSNVMODEL()
    if model_params:
        m = SCSNVMODEL(**model_params)

    model = m.model_construct()
    if hdf5_file and hdf5_fliename is not None:
        model.load_weights(hdf5_fliename)

    ## 用generator送数据
    if test_model == 1:
        testing_generator = DataGenerator(samples_test_data)
        if generator_params:
            testing_generator = DataGenerator(samples_test_data, **generator_params)

        loss, accuracy = model.evaluate_generator(generator=testing_generator, verbose=1)
        print('Accuracy: %f' % (accuracy * 100))
        result = model.predict_generator(generator=testing_generator, verbose=1)
        # print(result)
        return testing_generator.get_sendin_content()

    ## 直接送数据
    elif test_model == 2:
        loss, accuracy = model.evaluate(samples_test_data[0], samples_test_data[1])
        print('Accuracy: %f'%(accuracy * 100))
        result = model.predict(samples_test_data[0])
        return result


