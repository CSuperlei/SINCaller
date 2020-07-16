class TREGION:
    def __init__(self, sample_name):
        self.sample_name = sample_name

    def readfile(self, filename):
        region_file = open(filename)
        if region_file is not None:
            return region_file
        else:
            print('region file is empty')


    def region_info(self, region_file):
        ls = []
        lines = region_file.readlines()
        for line in lines:
            if line == '':
                continue
            line = line.strip()
            item = line.split(",")
            # print(item)
            sample = self.sample_name
            chr = item[0]
            left = item[1]
            right = item[2]
            # item = list(item)
            ls.append([sample, chr, left, right])

        return ls



