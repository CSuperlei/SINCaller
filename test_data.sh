#!/bin/bash
dir=/mnt/sdc/cailei/bio_project/scSNVIndel_data/test_data/test_data_list.txt

cat $dir | while read id
do
	python main.py  -ld /mnt/sdc/cailei/bio_project/scSNVIndel_data/test_data/${id}_wu_test_data.npy -sr /mnt/sdc/cailei/bio_project/scSNVIndel_data/predict_result_30/${id}_wu_predict_data.npy -g 1 -lo 3 -m 2 -tm 3
	python main.py  -ld /mnt/sdc/cailei/bio_project/scSNVIndel_data/predict_result_30/${id}_wu_predict_data.npy -fa /mnt/sdc/cailei/bio_common/human_ref/hg38.fa.fai -ov /mnt/sdc/cailei/bio_project/scSNVIndel_data/scSNVIndel_vcf_30/${id}_wu_predict.vcf -g 1 -lo 3 -m 5
done
