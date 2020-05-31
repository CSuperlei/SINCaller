from .scSNV_model import SCSNVMODEL


def testing(samples_test_data, samples_test_label, model_params=None, hdf5_file=True, hdf5_fliename='/home/cailei/bio_project/nbCNV/train_log/model_checkpoint/model_06.hdf5'):
    m = SCSNVMODEL()
    if model_params:
        m = SCSNVMODEL(**model_params)

    model = m.model_construct()
    if hdf5_file and hdf5_fliename is not None:
        model.load_weights(hdf5_fliename)

    loss, accuracy = model.evaluate(samples_test_data, samples_test_label, verbose=1)
    result = model.predict(samples_test_data[0])
    print('Accuracy: %f' % (accuracy * 100))
    print(result)