#!/bin/bash
dir=/mnt/sdc/cailei/bio_project/scSNVIndel_data/test_data/test_data_list.txt

cat $dir | while read id
do
	nohup python main.py -b /mnt/sdc/cailei/bio_project/SRP044380/bam_wu/${id}.sort.bam -fa /mnt/sdc/cailei/bio_common/human_ref/hg38.fa -r ./test_region.txt -sn ${id} -d /mnt/sdc/cailei/bio_project/scSNVIndel_data/test_data/${id}_wu_test_data.npy -g 1 -lo 1 -m 3 > ../scSNVIndel_data/nohup_data.out &
done

