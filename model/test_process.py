from .scSNV_model import SCSNVMODEL
from .train_data_generator import DataGenerator
from .test_data_generator import TEST

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
        testing_generator = TEST(samples_test_data, test_type=1)
        test_data = testing_generator.data_generator()
        print(test_data)
        loss, accuracy = model.evaluate(test_data, batch_size=64)
        print('Accuracy: %f'%(accuracy * 100))
        # print(testing_generator.get_sendin())

    ## 随机生成数据送入网络
    elif test_model == 3:
        testing_generator = TEST(samples_test_data, test_type=2)
        test_data = testing_generator.data_generator()
        result = model.predict(test_data, batch_size=64)
        print(result)
        print(testing_generator.get_sendin())


