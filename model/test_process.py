from .scSNV_model import SCSNVMODEL
from .train_data_generator import DataGenerator
from .test_data_generator import TEST
import numpy as np

def testing(samples_test_data, test_model=1, model_params=None,  generator_params=None, hdf5_file=True, hdf5_fliename='/home/cailei/bio_project/nbCNV/train_log/model_checkpoint/model_30.hdf5'):
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

        print(testing_generator)
        evaluate = model.evaluate_generator(generator=testing_generator, verbose=1)
        print(evaluate)
        print(testing_generator.get_sendin_content())

    ## 直接送有标签的测试数据进行评价
    elif test_model == 2:
        testing_generator = TEST(samples_test_data, test_type=1)
        X, Y = testing_generator.data_generator()
        evalue = model.evaluate(x=X, y=Y, batch_size=64)
        print(evalue)
        result = model.predict(x=X, batch_size=64)
        print(result)
        # re = np.array(result)
        # print(re.argmax())
        print(testing_generator.get_sendin())

    ## 直接送没有便签的测试数据进行评价
    elif test_model == 3:
        testing_generator = TEST(samples_test_data, test_type=2)
        test_data = testing_generator.data_generator()
        result = model.predict(test_data, batch_size=64)
        print(result)
        re_base = np.argmax(result[0], axis=1)
        print('re_base', re_base)
        re_indel = np.argmax(result[1], axis=1)
        print('re_indel', re_indel)
        re_genotype = np.argmax(result[2], axis=1)
        print('re_genotype', re_genotype)

        print(testing_generator.get_sendin())


