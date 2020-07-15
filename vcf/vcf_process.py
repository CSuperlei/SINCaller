from pysam import VariantFile
from textwrap import dedent

class VCF:
    def readfile(self, filename):
        bcf_in = VariantFile(filename, 'r')
        if None == bcf_in:
            print('vcf file is empty')
   
        return bcf_in

    def varient_info(self, vcf_file):
        ls = []
        for rec in vcf_file.fetch():
            sample = ''
            label = ()
            for key, value in rec.samples.items():
                sample = key
                label = value['GT']
                break
            s_c_p = sample + '_' + rec.chrom + '_' + str(rec.pos)
            # print(rec.pos, rec.ref, rec.alts[0])
            ref = list(rec.ref)
            alts = list(rec.alts[0])
            # print(rec.pos, ref[-1], alts[-1])
            ls.append((s_c_p, (ref[-1], alts[-1]), label))
        return ls

    def gernate_vcf_title(self, title, fastai, sample_name, output_file):
        if output_file is not None:
            output_file = open(output_file, 'w')
            def output(string):
                print(string, file=output_file)
        title = dedent(title)
        output(title)
        if fastai is not None:
            with open(fastai, "r") as fai_fp:
                for row in fai_fp:
                    columns = row.strip().split("\t")
                    contig_name, contig_size = columns[0], columns[1]
                    output("##contig=<ID=%s,length=%s>" % (contig_name, contig_size)) ## 染色体名称和长度
        output('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % (sample_name))

    def generate_vcf_content(self, output_file, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, VALUE):
        if output_file is not None:
            with open(output_file, "a") as output_file:
                def output(string):
                    print(string, file=output_file)

                print(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, VALUE)
                output(CHROM + '\t' + POS + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' + FILTER + '\t' + INFO + '\t' + FORMAT + '\t' + VALUE)







