
class TREGION:
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
            item = line.split(",")
            print(item)
            sample = item[0]
            chr = item[1]
            left = item[2]
            right = item[3]
            # item = list(item)
            ls.append([sample, chr, left, right])

        return ls



