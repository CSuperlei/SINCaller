import os
import numpy as np

class DATACOMBINE:
    def __init__(self, path, target):
        self.path = path
        self.target = target

    def data_combine(self):
        files = os.listdir(self.path) #得到文件夹下的所有文件名称
        re = []
        for file in files: #遍历文件夹
            # print(file)
            if not os.path.isdir(file): #判断是否是文件夹，不是文件夹才打开
                data = np.load(self.path+"/"+file, allow_pickle=True)
                data = data.tolist()
                re += data

        print(re)
        np.save(self.target, re)


def main():
    # path = "f:\\Research\\nbCNV\\train_data\\"  # 文件夹目录
    # target = "f:\\Research\\nbCNV\\train_summary_data\\training_data.npy"
    path = '/home/cailei/bio_project/nbCNV/train_data'
    target = '/home/cailei/bio_project/nbCNV/train_summary_data/training_data.npy'
    dc = DATACOMBINE(path, target)
    dc.data_combine()


if __name__ == '__main__':
    main()