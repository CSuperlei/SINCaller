import numpy as np

def data_bug(dir):
    data_sample = np.load(dir)

    for i in data_sample:
        if '' in set(i[1][1]):
            print(i)
            re = ['d' if item == '' else item for item in i[1][1]]
            print(re)


if __name__ == '__main__':
    dir = 'f:\\Research\\Bio_Project\\nbCNV\\train_data\\train_data.npy'
    data_bug(dir)



