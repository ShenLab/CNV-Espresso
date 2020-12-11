from __future__ import division
import pysam
import vcf
import os,sys
import gzip
import pdb

def chrStr2Int(chr):
    if type(chr) == str:
        chr = chr.replace('chr','')
    else:
        pass
    if chr == 'X' or chr == 'x':
        return 23
    elif chr == 'Y' or chr == 'y':
        return 24
    else:
        return int(chr)

def open_capture_region(file):
    fp = open(file,'r')

    capture_region_dict ={}
    for one_line in fp:
        one_line = one_line.strip()
        one_line_list = one_line.split('\t')

        exon_chr = chrStr2Int(one_line_list[0])
        exon_start = int(one_line_list[1])
        exon_stop = int(one_line_list[2])

        if capture_region_dict.has_key(exon_chr):
            capture_region_dict[exon_chr].append([exon_chr,exon_start,exon_stop])
        else:
            capture_region_dict[exon_chr] = [[exon_chr,exon_start,exon_stop]]

    return capture_region_dict

def overlap_region(region1_start,region1_stop,region2_start,region2_stop):
    if region1_start == region2_start and region1_stop == region2_stop:
        case = 'Case 6'
        overlap_region_start = region1_start
        overlap_region_stop = region1_stop
    elif region1_start < region2_start and region1_stop >= region2_start and region1_stop <= region2_stop:
        case = 'Case 1'
        overlap_region_start = region2_start
        overlap_region_stop = region1_stop
    elif region1_start >= region2_start and region1_start <= region2_stop and region1_stop > region2_stop:
        case = 'Case 2'
        overlap_region_start = region1_start
        overlap_region_stop = region2_stop
    elif region1_start >= region2_start and region1_stop <= region2_stop:
        case = 'Case 3'
        overlap_region_start = region1_start
        overlap_region_stop = region1_stop
    elif region2_start>= region1_start and region2_stop <= region1_stop:
        case = 'Case 4'
        overlap_region_start = region2_start
        overlap_region_stop = region2_stop
    else:
        case = 'Case 5 or 7, No overlap.'
        overlap_region_start = None
        overlap_region_stop = None

    return ([overlap_region_start,overlap_region_stop,case])

def find_overlap_exon(overlap_region_list,exome_region_dict,output_unique,overlap_threshold):
    overlap_exon_list = []
    overlap_exon_list_unique = []

    for overlap_region_reader in overlap_region_list:
        cnp_name = overlap_region_reader[0]
        cnp_chr = chrStr2Int(overlap_region_reader[1])
        cnp_start = int(overlap_region_reader[2])
        cnp_stop = int(overlap_region_reader[3])
        snp_pos = int(overlap_region_reader[4])
        R2 = float(overlap_region_reader[5])
        population = overlap_region_reader[6]

        try:
            for exon in exome_region_dict[cnp_chr]:
                exon_chr = chrStr2Int(exon[0])
                exon_start = int(exon[1]) + 1 #for bed file
                exon_stop = int(exon[2]) + 1 #for bed file
                exon_size = exon_stop - exon_start + 1

                if cnp_chr == exon_chr:
                    result_start = None
                    result_stop = None
                    result_case = None

                    result_overlap_region = overlap_region(cnp_start,cnp_stop,exon_start,exon_stop)
                    result_start = result_overlap_region[0]
                    result_stop = result_overlap_region[1]
                    result_case = result_overlap_region[2]

                    if result_start != None and result_stop != None:
                        overlap_ratio = (result_stop - result_start + 1)/exon_size
                        if overlap_ratio >= overlap_threshold:
                            overlap_exon_list.append([cnp_name,cnp_chr,cnp_start,cnp_stop,snp_pos,R2,population,exon_chr,exon_start,exon_stop])
                            print cnp_name,cnp_chr,cnp_start,cnp_stop,snp_pos,R2,population,exon_chr,exon_start,exon_stop
                            #store the unique data
                            if [exon_chr,exon_start,exon_stop] not in overlap_exon_list_unique:
                                overlap_exon_list_unique.append([exon_chr,exon_start,exon_stop])
        except KeyError:
            print 'dict Key Error!, No chromosome %s in exome_region_dict'%cnp_chr
        except:
            print 'other errors!'
            return

    if output_unique == True:
        '''Return uniqued [exon_chr,exon_start,exon_stop]'''
        return overlap_exon_list_unique
    elif output_unique == False:
        '''Return  [sample_name,status,exon_chr,exon_start,exon_stop]'''
        return overlap_exon_list
    else:
        print 'Parement output_unique error.'
        return

def output_to_file(results,output_file):
    fp = open(output_file,'w')
    if type(results) == list:
        for result_line in results:
            for each_one in result_line:
                fp.write(str(each_one))
                fp.write('\t')
            fp.write('\n')

    elif type(results) == dict:
        for key in results.keys():
            for result_line in results[key]:
                for each_one in result_line:
                    fp.write(str(each_one))
                    fp.write('\t')
                fp.write('\n')
    else:
        print 'unsupport results type!'
    print 'The variable has already output to %s\n\n'%output_file
    fp.close()


tagSNP_file='tagSNP_hg19.txt'
vcf_file_gz='NA12878.vcf.gz'
TargetRegion_file = '20130108.exome.targets_noChr.bed'
output_file='tagSNP_results.txt'
input_population = 'CEU'

TargetRegion_dict = open_capture_region(TargetRegion_file)

vcf_reader = vcf.Reader(filename=vcf_file_gz)
vcf_file_name = os.path.basename(vcf_file_gz)
vcf_file_name = vcf_file_name.split('.',1)[0]
tagsnp_reader=open(tagSNP_file,'r')


i=0
CNP_list = []
for row in tagsnp_reader:
    i=i+1
    if i != 1:
        row=row.strip()
        element= row.split('\t')
        cnp_name= element[0]
        chr= element[1].replace("chr","")
        cnp_start=element[2]
        cnp_stop=element[3]
        snp_pos=int(element[4])
        R2 = float(element[5])
        population = element[6]
 
        if input_population == population:
            for vcf_snp in vcf_reader.fetch(chr,snp_pos-1,snp_pos):
                if len(str(vcf_snp)) !=0:
                    CNP_list.append([cnp_name,chr,cnp_start,cnp_stop,snp_pos,R2,population])

print 'CNP_list number:',len(CNP_list)
CNP_exon_list = find_overlap_exon(CNP_list,TargetRegion_dict,output_unique=False,overlap_threshold=0.5)
print 'CNP_exon number:',len(CNP_exon_list)
output_to_file(CNP_exon_list,output_file)

tagsnp_reader.close()
