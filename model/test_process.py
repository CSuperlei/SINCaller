from .scSNV_model import SCSNVMODEL


def testing(samples_test_data, model_params=None, hdf5_file=True, hdf5_fliename='/home/cailei/bio_project/nbCNV/train_log/model_checkpoint/model_08.hdf5'):
    m = SCSNVMODEL()
    if model_params:
        m = SCSNVMODEL(**model_params)

    model = m.model_construct()
    if hdf5_file and hdf5_fliename is not None:
        model.load_weights(hdf5_fliename)

    loss, accuracy = model.evaluate(samples_test_data[0], samples_test_data[1], verbose=1)
    print('Accuracy: %f' % (accuracy * 100))