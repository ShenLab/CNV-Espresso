from __future__ import division
import vcf, os
import csv
import re
import pdb
import datetime
import numpy as np
import pysam
import sys
import bz2

def main():
    pass
# input/output/parse files
def parse_pedigreeFile(file_name):
    unaffected_parents_list = []
    try:
        fp = open(file_name,'r')
        row_num = 0
        ped_indiv_dict = {}
        ped_family_dict = {}

        for reader in fp:
            row_num += 1
            reader = reader.strip()
            reader = reader.split('\t')
            FamilyID = reader[0]
            IndividualID = reader[1]
            PaternalID = reader[2]
            MaternalID = reader[3]
            Sex = reader[4]
            Phenotype = reader[5]
            if len(reader) == 7:
                Role = reader[6]
                ped_indiv_dict[IndividualID] = [FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role]
                if FamilyID in ped_family_dict:
                    ped_family_dict[FamilyID].append([FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role])
                else:
                    ped_family_dict[FamilyID] = [[FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role]]
            elif len(reader) >= 8:
                Role = reader[6]
                Ancestry = reader[7]
                ped_indiv_dict[IndividualID] = [FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role, Ancestry]
                if FamilyID in ped_family_dict:
                    ped_family_dict[FamilyID].append([FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role, Ancestry])
                else:
                    ped_family_dict[FamilyID] = [[FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role, Ancestry]]
            else:
                ped_indiv_dict[IndividualID] = [FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype]
                if FamilyID in ped_family_dict:
                    ped_family_dict[FamilyID].append([FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype])
                else:
                    ped_family_dict[FamilyID] = [[FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype]]
        
        time_stamp = datetime.datetime.now()
        print(time_stamp.strftime('%Y.%m.%d-%H:%M:%S'),"\tPedigree file %s imported."%file_name)
        return [ped_indiv_dict, ped_family_dict]
    except:
        print('Error df.parse_pedigreeFile function when openning %s'%file_name)
        return [None, None]

def parseXhmmXcnvFile(xcnv_file):
    header_list = getHeader(xcnv_file)
    sample_id_col = header_list.index(fetch_colName(['SAMPLE','sample_ID'],header_list))
    cnv_type_col = header_list.index(fetch_colName(['cnv_type','type','CNV'],header_list))
    cnv_interval_col = header_list.index(fetch_colName(['INTERVAL'],header_list))
    sq_col = header_list.index(fetch_colName(['Q_SOME'],header_list))

    fp = open(xcnv_file,'r')
    pattern = re.compile(r'(\S+):(\d+)-(\d+)')
    cnv_num = 0
    cnv_list=[]
    for row in fp:
        cnv_num += 1
        if cnv_num == 1:
            pass
        else:
            row = row.strip()
            row = row.split('\t')
            sample_id = row[sample_id_col]
            sample_id = sample_id.replace('NG-','')
            sample_id = sample_id.replace('-001-001','')
            sample_id = sample_id.replace('-002-002','')
            sample_id = sample_id.replace('-003-003','')
            sample_id = sample_id.replace('-004-004','')
            sample_id = sample_id.replace('-005-005','')
            sample_id = sample_id.replace('-006-006','')
            sample_id = sample_id.replace('-007-007','')
            sample_id = sample_id.replace('-008-008','')
            sample_id = sample_id.replace('-009-009','')
            sample_id = sample_id.replace('-010-010','')
            sample_id = sample_id.replace('-011-011','')
            sample_id = sample_id.replace('-012-012','')
            sample_id = sample_id.replace('-013-013','')
            sample_id = sample_id.replace('-014-014','')
            sample_id = sample_id.replace('-015-015','')
            sample_id = sample_id.replace('-016-016','')

            cnv_type = row[cnv_type_col]
            location = row[cnv_interval_col]
            loc_list = pattern.match(location)
            cnv_chr = loc_list.group(1)
            cnv_start = int(loc_list.group(2))
            cnv_stop = int(loc_list.group(3))
            if row[sq_col] != 'NA':
                cnv_SQ = float(row[sq_col])
            else:
                cnv_SQ = -1
            cnv_list.append([sample_id, cnv_chr, cnv_start, cnv_stop, cnv_type, cnv_SQ])
    fp.close()
    return cnv_list

def parseInterval(cnv_interval):
    pattern = re.compile(r'(\S+):(\d+)-(\d+)')
    try:
        loc_list = pattern.match(cnv_interval)
        cnv_chr = loc_list.group(1)
        cnv_start = int(loc_list.group(2))
        cnv_stop = int(loc_list.group(3))
    except:
        print("[ERROR] interval does not match 'chr:start-stop' pattern.")
        cnv_chr = None
        cnv_start = None
        cnv_stop = None
    return [cnv_chr, cnv_start, cnv_stop]

def parsing_ExAC_CNV_file(ExAC_cnv_file):
    fp = open(ExAC_cnv_file,'r')
    ExAC_cnv_num = 0
    ExAC_cnv_list=[]

    for row in fp:
        row = row.strip()
        row = row.split('\t')
        
        chromesome = row[0].replace("chr","")
        start = int(row[1])
        stop = int(row[2])
        status = row[3]
        length = stop - start + 1
        ExAC_cnv_list.append([chromesome,start,stop,status])
        ExAC_cnv_num += 1

    fp.close()
    return ExAC_cnv_list

def import_stored_results(input_file):
    if(os.path.isfile(input_file) is False):
        path = os.getcwd()
        input_file = path + input_file

    time_stamp = datetime.datetime.now()
    print(time_stamp.strftime('%Y.%m.%d-%H:%M:%S'),"Importing already dictionary file to Memory....")

    result_file_dict = {}
    fp = open(input_file,'r')
    for reader in fp:
        reader = reader.strip()
        reader_list = reader.split('\t')
        key = reader_list[0]
        if key in result_file_dict:
            result_file_dict[key].append(reader_list)
        else:
            result_file_dict[key] = [reader_list]
    fp.close()
    print(input_file,' imported!')
    return result_file_dict

def parse_pseq_denovo_file(cnv_file):
    fp = open(cnv_file,'r')
    pattern = re.compile(r'(\S+):(\d+)..(\d+)')
    pattern_cnv_type = re.compile(r'<(\S+)>')
    cnv_num = 0
    cnv_list=[]
    for row in fp:
        if cnv_num == 0:
            cnv_num += 1
            cnv_list.append(['status','sample_name','cnv_chr','cnv_start','cnv_stop','cnv_length','cnv_type','num_target'])
        else:
            row = row.strip()
            row = row.split('\t')
            status = row[0]
            sample_name = row[1]
            location = row[2]
            loc_list = pattern.match(location)
            cnv_chr = loc_list.group(1).replace('chr','')
            cnv_start = int(loc_list.group(2))
            cnv_stop = int(loc_list.group(3))
            cnv_length = cnv_stop - cnv_start + 1
            cnv_type = pattern_cnv_type.match(row[3]).group(1)
            num_target = int(row[4])
            cnv_list.append([status,sample_name,cnv_chr,cnv_start,cnv_stop,cnv_length,cnv_type,num_target])
            # print(cnv_num,status,sample_name,cnv_chr,cnv_start,cnv_stop,cnv_type,num_target)
            cnv_num += 1
    fp.close()
    return cnv_list

def output_to_file(results,file_name):
    path = os.path.dirname(file_name)
    if not os.path.exists(path) and path!='':
        os.mkdir(path)

    fp = open(file_name,'w')
    if type(results) == list:
        for result_line in results:
            result_line = [str(x) for x in result_line]
            fp.write("\t".join(result_line) + os.linesep)            

    elif type(results) == dict:
        for key in results.keys():
            for result_line in results[str(key)]:
                result_line = [str(x) for x in result_line]
                fp.write("\t".join(result_line) + os.linesep)

    else:
        print('[ERROR]: Unsupport results type!')
        return
    
    print('[INFO]: File outputs to %s'%file_name)
    time_stamp = datetime.datetime.now()
    print(time_stamp.strftime('[%Y.%m.%d-%H:%M:%S]'))
    fp.close()

def vcfToDict(vcf_file):
    vcf_reader = vcf.Reader(filename = vcf_file)
    sample_list = vcf_reader.samples
    variant_dict = {}
    num = 0

    time_stamp = datetime.datetime.now()
    print(time_stamp.strftime('%Y.%m.%d-%H:%M:%S'),"\tParsing VCF file to dictionary....")
    for record in vcf_reader: 
        num += 1
        if num % 1000 == 1:
            print(num, "variant_dict:",len(variant_dict)) 
            print(record)
    
        for sample in sample_list:
            DSCVR = record.genotype(sample)['DSCVR']
            if DSCVR == "Y":
                # print(sample,'---------------------')
                # print(record)
                # print(record.INFO)
                # print(record.genotype(sample))
                cnv_chr = record.CHROM
                cnv_start = record.INFO['TPOS']
                cnv_stop = record.INFO['END']
           
                if record.genotype(sample)['GT'] == '.':
                    cnv_type = None
                else:
                    GT = int(record.genotype(sample)['GT'])
                    # print(GT)
                    if GT == 0:
                        cnv_type = record.REF
                    else:
                        cnv_type = record.ALT[GT-1]
                        AF = record.INFO['AF'][GT-1]
                        NUMT = record.INFO['NUMT']
                        SQ = record.genotype(sample)['SQ'][GT-1]
                        EQ = record.genotype(sample)['EQ'][GT-1]
                        NQ = record.genotype(sample)['NQ'][GT-1]
                        NDQ = record.genotype(sample)['NDQ']

                        if sample in variant_dict:
                            variant_dict[sample].append([sample, cnv_chr, cnv_start, cnv_stop, cnv_type, SQ, EQ, NQ, NDQ, AF, NUMT])
                        else:
                            variant_dict[sample] = [[sample, cnv_chr, cnv_start, cnv_stop, cnv_type, SQ, EQ, NQ, NDQ, AF, NUMT]]
    return variant_dict

def fileToDict(file_name, start_row=0, key_col=0, tab_sep_only = False):
    result_dict = {}
    fp = open(file_name)
    num = 0
    for row in fp:
        num += 1
        if num < start_row:
            continue
        row = row.strip()

        if tab_sep_only == True:
            row = re.split('\t| ',row)
        else:
            row = re.split(r'[;,\s]\s*', row)

        row = [s.replace('"','') for s in row]
        if row[key_col] in result_dict:
            result_dict[row[key_col]].append(row)         
        else:
            result_dict[row[key_col]] = [row]
    print("Imported %s."%file_name)
    return result_dict

def xhmmRDToDict(RD_file):
    num = 0
    rd_dict = {}
    print('Importing %s'%RD_file)
    with open(RD_file, "rb") as f:
        for reader in f:
            num += 1
            reader = reader.strip()
            reader = reader.split()
            sampleID = reader[0].decode("utf-8")
            sampleID = sampleID.replace('NG-','')
            sampleID = sampleID.replace('-001-001','')
            sampleID = sampleID.replace('-002-002','')
            sampleID = sampleID.replace('-003-003','')
            sampleID = sampleID.replace('-004-004','')
            sampleID = sampleID.replace('-005-005','')
            sampleID = sampleID.replace('-006-006','')
            sampleID = sampleID.replace('-007-007','')
            sampleID = sampleID.replace('-008-008','')
            sampleID = sampleID.replace('-009-009','')
            sampleID = sampleID.replace('-010-010','')
            sampleID = sampleID.replace('-011-011','')
            sampleID = sampleID.replace('-012-012','')
            sampleID = sampleID.replace('-013-013','')
            sampleID = sampleID.replace('-014-014','')
            sampleID = sampleID.replace('-015-015','')
            sampleID = sampleID.replace('-016-016','')
            rd_dict[sampleID] = reader[1:]
            if num%1000 == 0:
                print('Importing...',num,sampleID)
    return(rd_dict)

def xhmmRDToDict_bySampleID(RD_file, Child_ID, Father_ID, Mother_ID):
    num = 0
    RD_dict = {}
    if Father_ID != None and Mother_ID != None:
        command = 'cat %s | grep "Matrix\|%s\|%s\|%s"'%(RD_file, Child_ID, Father_ID, Mother_ID)
    else:
        command = 'cat %s | grep "Matrix\|%s"'%(RD_file, Child_ID)

    print(command)
    stream = os.popen(command)
    RD_list = stream.readlines()
    for RD_reader in RD_list:
        RD_reader = RD_reader.strip().split()
        # sampleID = RD_reader[0].decode("utf-8")
        sampleID = str(RD_reader[0])
        sampleID = sampleID.replace('NG-','')
        sampleID = sampleID.replace('-001-001','')
        sampleID = sampleID.replace('-002-002','')
        sampleID = sampleID.replace('-003-003','')
        sampleID = sampleID.replace('-004-004','')
        sampleID = sampleID.replace('-005-005','')
        sampleID = sampleID.replace('-006-006','')
        sampleID = sampleID.replace('-007-007','')
        sampleID = sampleID.replace('-008-008','')
        sampleID = sampleID.replace('-009-009','')
        sampleID = sampleID.replace('-010-010','')
        sampleID = sampleID.replace('-011-011','')
        sampleID = sampleID.replace('-012-012','')
        sampleID = sampleID.replace('-013-013','')
        sampleID = sampleID.replace('-014-014','')
        sampleID = sampleID.replace('-015-015','')
        sampleID = sampleID.replace('-016-016','')
        RD_dict[sampleID] = RD_reader[1:]
    return(RD_dict)

def fileToList(file_name):
    result_list = [] 
    extension = os.path.splitext(file_name)[-1][1:]
    if extension == 'bz2':
        fp = bz2.open(file_name, "rt")
    else:
        fp = open(file_name)
    for row in fp:
        row = row.strip()
        row = row.replace("\"","")
        #row = re.split(r'[;,\s]\s*', row)
        row = row.split('\t')
        result_list.append(row)
    fp.close()  
    return result_list 

def fileToList_tab(file_name):
    result_list = [] 
    fp = open(file_name)
    for row in fp:
        row = row.strip()
        row = row.replace("\"","")
        row = re.split('\t| ',row)
        #row = re.split(r'[;,\s]\s*', row)
        result_list.append(row)
    fp.close()  
    return result_list 

def printList(list_name):
    num = 0
    for reader in list_name:
        num += 1
        print(num, reader)

# fetch
def fetchVCF_BAF(vcf_file, sample_id, cnv_chr, cnv_start, cnv_stop):
    vcf_reader = vcf.Reader(filename=vcf_file)
    samples_list = vcf_reader.samples
    if sample_id not in samples_list:
        return ['NA', 'NA', 'NA', 'NA']
    else:
        records = vcf_reader.fetch(cnv_chr, cnv_start, cnv_stop + 1)
        baf_list = []
        snp_num, het_num = 0, 0    
        for record in records:
            Zygosity, baf, DP = None, None, None

            if record.is_snp and record.genotype(str(sample_id)).is_variant: 
                #     print(snp_num, record, record.genotype(str(sample_id)), record.is_snp, record.genotype(str(sample_id)).is_variant)
                #     pdb.set_trace()
                snp_num += 1
                if record.genotype(str(sample_id)).is_het:
                    het_num += 1
                    Zygosity = 'Het'
                else:
                    Zygosity = 'Hom'

                sample_snv = record.genotype(str(sample_id)).data
                DP = sample_snv.DP
                ref_allele_reads = sample_snv.AD[0]
                alt_allele_reads = sum(x for x in sample_snv.AD[1: len(sample_snv.AD)])
                try:
                    baf = round(alt_allele_reads/(ref_allele_reads+alt_allele_reads),2)
                except Exception as e:               
                    continue
                print(record, sample_snv, Zygosity, baf, DP)
                baf_list.append(Zygosity+'|'+str(baf)+'|'+str(DP))     

        baf_info = ",".join(baf_list)
        if baf_info == '':
            baf_info = '-'

        if snp_num != 0:
            het_ratio = round(het_num/snp_num,2)
        else:
            het_ratio = '-'
        return [snp_num, baf_info, het_num, het_ratio]

def fetchVCF_SQ_NQ(vcf_file, sampleID, cnv_chr, cnv_start, cnv_stop, cnv_type):
    vcf_reader = vcf.Reader(filename=vcf_file)
    samples_list = vcf_reader.samples
    if sampleID not in samples_list:
        return [None, None]
    else:
        SQ, NQ = None, None
        try:
            for record in vcf_reader.fetch(cnv_chr, cnv_start, cnv_stop + 1):
                if record.INFO['TPOS'] == cnv_start and record.INFO['END'] == cnv_stop:
                    if str(cnv_type).upper() in str(record.ALT[0]):
                        GT = 0
                    elif str(cnv_type).upper() in str(record.ALT[1]):
                        GT = 1
                    else:
                        pdb.set_trace()
                    try:
                        SQ = record.genotype(sampleID)['SQ'][GT]
                        NQ = record.genotype(sampleID)['NQ'][GT]
                    except:
                        SQ = record.genotype(sampleID)['SQ']
                        NQ = record.genotype(sampleID)['NQ']
                    break
        except:
            pdb.set_trace()
        return [SQ, NQ]

def fetchCNVs(family_ID, sampleID, vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file,'r')) 
    cnv_list = []
    for record in vcf_reader:
        if record.genotype(sampleID)['DSCVR'] == 'Y':
            cnv_chr = record.CHROM
            cnv_start = record.INFO['TPOS']
            cnv_stop = record.INFO['END']
            if record.genotype(sampleID)['GT'] == '1':
                cnv_type = 'DEL'
            if record.genotype(sampleID)['GT'] == '2':
                cnv_type = 'DUP'
            SQ = max(record.genotype(sampleID)['SQ'])
            EQ = max(record.genotype(sampleID)['EQ'])
            NQ = max(record.genotype(sampleID)['NQ'])
            cnv_list.append([family_ID, sampleID, cnv_chr, cnv_start, cnv_stop, cnv_type, SQ, EQ, NQ])
    return cnv_list   

def fetch_cnv_vcf(vcf_file, sample_id, cnv_chr, cnv_start, cnv_stop):
    vcf_reader = vcf.Reader(filename=vcf_file)
    records = vcf_reader.fetch(cnv_chr, cnv_start, cnv_stop+1)
    mean_rd = '-'
    mean_orig_rd = '-'
    for record in records:
        if record.genotype(sample_id)['DSCVR'] == 'Y':
            # print(record.INFO, record.genotype(sample_id))
            mean_rd = str(record.genotype(sample_id).data.RD)
            mean_orig_rd = str(record.genotype(sample_id).data.ORD)
    return [mean_rd, mean_orig_rd]

def fetchAllelFrequencyinVCF(vcf_file, chr, start ,stop):
    vcf_reader = vcf.Reader(filename = vcf_file)
    start = int(start)
    stop = int(stop)
    for record in vcf_reader.fetch(chr,start,stop):
        if record.INFO['TPOS'] == start and record.INFO['TEND'] == stop:
            af_max = max(record.INFO['AF'])
    return af_max

def fetchClammsSpecialRegions(clamms_special_regions_list, cnv_chr, cnv_start, cnv_stop):
    result = []
    for reader in clamms_special_regions_list:
        clamms_special_region_chr = reader[0]
        clamms_special_region_start = int(reader[1])
        clamms_special_region_stop = int(reader[2])
        cnv_start = int(cnv_start)
        cnv_stop = int(cnv_stop)
        if regionOverlap_occupyFirst(cnv_chr, cnv_start, cnv_stop, clamms_special_region_chr, 
            clamms_special_region_start, clamms_special_region_stop, 0.5) == True:
            result = clamms_special_region_chr+':'+str(clamms_special_region_start)+'-'+str(clamms_special_region_stop)
            break

    return result

def fetchSegmentalDuplicationRegions(SegmentalDuplication_list, cnv_chr, cnv_start, cnv_stop):
    result = []
    for reader in SegmentalDuplication_list:
        sd_chr = reader[0].replace("chr","")
        sd_start = int(reader[1])
        sd_stop = int(reader[2])
        cnv_start = int(cnv_start)
        cnv_stop = int(cnv_stop)
        if regionOverlap_occupyFirst(cnv_chr, cnv_start, cnv_stop, sd_chr, sd_start, sd_stop,  0.75) == True:
            result = 'SD_'+sd_chr+':'+str(sd_start)+'-'+str(sd_stop)
            break
    return result

def fetchSegmentalDuplicationRegions_Dict(SegmentalDuplication_dict, cnv_chr, cnv_start, cnv_stop):
    result = []
    cnv_chr = cnv_chr.replace("chr","")
    for reader in SegmentalDuplication_dict["chr"+cnv_chr]:
        sd_chr = reader[0].replace("chr","")
        sd_start = int(reader[1])
        sd_stop = int(reader[2])
        cnv_start = int(cnv_start)
        cnv_stop = int(cnv_stop)
        if regionOverlap_occupyFirst(cnv_chr, cnv_start, cnv_stop, sd_chr, sd_start, sd_stop,  0.75) == True:
            result = 'SD_'+sd_chr+':'+str(sd_start)+'-'+str(sd_stop)
            break
    return result

def fetchExAC_CNV(ExAC_list, cnv_chr, cnv_start, cnv_stop):
    result = []
    for reader in ExAC_list:
        try:
            region_chr = reader[0].replace("chr","")
            region_start = int(reader[1])
            region_stop = int(reader[2])
            cnv_start = int(cnv_start)
            cnv_stop = int(cnv_stop)
            cnv_type = None
            region_type = None            
            if reciprocalOverlap(cnv_chr, cnv_start, cnv_stop, cnv_type, region_chr, region_start, region_stop, region_type, 0.5) == True:
                result = reader
                break
        except:
            pass

    return result

def fetch_gnomAD_CNV(gnomAD_list, cnv_chr, cnv_start, cnv_stop, cnv_type):
    result_list = []
    for reader in gnomAD_list:
        try:
            gnomAD_chr = reader[0].replace("chr","")
            gnomAD_start = int(reader[1])
            gnomAD_stop = int(reader[2])
            gnomAD_type = reader[3]
            gnomAD_AF = reader[4]
            cnv_start = int(cnv_start)
            cnv_stop = int(cnv_stop)
            if gnomAD_type == cnv_type:
                if reciprocalOverlap(cnv_chr, cnv_start, cnv_stop, cnv_type, gnomAD_chr, gnomAD_start, gnomAD_stop, gnomAD_type, 0.5) == True:
                    result_list = [gnomAD_chr +':'+ str(gnomAD_start) +'-'+ str(gnomAD_stop) +'_'+ gnomAD_type +'_'+ gnomAD_AF, gnomAD_AF]
                    break
        except:
            pass
    return result_list

def fetch_spark_CNV(spark_list, cnv_chr, cnv_start, cnv_stop, cnv_type):
    result = []
    for reader in spark_list:
        try:
            spark_chr = reader[0].replace("chr","")
            spark_start = int(reader[1])
            spark_stop = int(reader[2])
            spark_type = reader[3]
            spark_sampleID = reader[4]
            
            cnv_start = int(cnv_start)
            cnv_stop = int(cnv_stop)
            if spark_type == cnv_type:
                if reciprocalOverlap(cnv_chr, cnv_start, cnv_stop, cnv_type, spark_chr, spark_start, spark_stop, spark_type, 0.5) == True:
                    result = spark_chr+':'+str(spark_start)+'-'+str(spark_stop)+'_'+spark_type+'_'+spark_sampleID
                    break
        except:
            pass

    return result

def fetchGene_scores(Gene_list, cnv_chr, cnv_start, cnv_stop, coordinate):
    gene_result = []
    pLI_result = []
    s_het_result = []
    Episcore_result = []
    num = 0
    for reader in Gene_list:
        num += 1
        if num == 1:
            pass
        else:
            region_chr = reader[0].replace("chr","")
            region_start = int(reader[1])
            region_stop = int(reader[2])
            gene_name = reader[3]

            pLI = reader[4]
            if pLI != '-':
                pLI = round(float(pLI),2)

            s_het = reader[5]
            if s_het != '-':
                s_het = round(float(s_het),2)

            Episcore = reader[6]
            if Episcore != '-':
                Episcore = round(float(reader[6]),2)

            cnv_start = int(cnv_start)
            cnv_stop = int(cnv_stop)

            if regionOverlap_occupyFirst(cnv_chr, cnv_start, cnv_stop, region_chr, region_start, region_stop, 1) == True:
                gene_info = gene_name+'_'+region_chr+':'+str(region_start)+'-'+str(region_stop)
                if coordinate == True:
                    gene_result.append(gene_info)
                else:
                    gene_result.append(gene_name)
                
                if pLI == '-':
                    pass
                elif pLI > 0.5:
                    pLI_result.append(gene_name+'='+str(pLI))

                if s_het == '-':
                    pass
                elif s_het > 0.09:
                    s_het_result.append(gene_name+'='+str(s_het))
                
                if Episcore == '-':
                    pass
                elif Episcore > 0.6:
                    Episcore_result.append(gene_name+'='+str(Episcore))

    if coordinate == False:
        gene_result = np.unique(gene_result)
        gene_result = ",".join(gene_result)

    pLI_result = np.unique(pLI_result)
    pLI_result = ",".join(pLI_result)
    s_het_result = np.unique(s_het_result)
    s_het_result = ",".join(s_het_result)
    Episcore_result = np.unique(Episcore_result)
    Episcore_result = ",".join(Episcore_result)

    return [gene_result, pLI_result, s_het_result, Episcore_result]

def fetchGeneScores_byGeneName(gene, Gene_list):
    num = 0
    for reader in Gene_list:
        num += 1
        if num == 1:
            pass
        else:
            pLI_result = ""
            s_het_result = ""
            Episcore_result = ""

            gene_name = reader[3]

            pLI = reader[4]
            if pLI != '-':
                pLI = round(float(pLI),2)

            s_het = reader[5]
            if s_het != '-':
                s_het = round(float(s_het),2)

            Episcore = reader[6]
            if Episcore != '-':
                Episcore = round(float(reader[6]),2)

            if str.lower(gene) == str.lower(gene_name):  
                # pdb.set_trace()              
                if pLI == '-':
                    pass
                elif pLI > 0.5:
                    pLI_result = str(pLI)

                if s_het == '-':
                    pass
                elif s_het > 0.09:
                    s_het_result = str(s_het)
                
                if Episcore == '-':
                    pass
                elif Episcore > 0.6:
                    Episcore_result = str(Episcore)
                break

    return [gene, pLI_result, s_het_result, Episcore_result]

def getParentsID(sampleID, ped_file):
    [ped_indiv_dict, ped_family_dict] = parse_pedigreeFile(ped_file)
    paternalID = ped_indiv_dict[sampleID][2]
    maternalID = ped_indiv_dict[sampleID][3]
    return([paternalID, maternalID])

def getBatchInfo(sampleID, batch_file):
    batch_dict = fileToDict(batch_file)
    return(batch_dict[sampleID][0][2])

def fetch_colName(keyWord_list, colName_list):
    for keyWord in keyWord_list:
        if keyWord in colName_list:
            return keyWord
        else:
            keyWord = None
    return keyWord

def getHeader(input_file):
    header_list = []
    fp = open(input_file)
    for reader in fp:
        reader = reader.strip()
        reader = re.split('\t| ',reader)
        header_list = reader
        break
    fp.close()
    return header_list

def fetch_cytoband(cytoband_list, cnv_chr, cnv_start, cnv_stop):
    result = []
    for reader in cytoband_list:
        try:
            cytoband_chr = reader[0].replace("chr","")
            cytoband_start = int(reader[1])
            cytoband_stop = int(reader[2])
            cytoband_info = reader[3]
            
            cnv_start = int(cnv_start)
            cnv_stop = int(cnv_stop)

            if regionOverlap_occupyFirst(cnv_chr, cnv_start, cnv_stop, cytoband_chr, cytoband_start, cytoband_stop, 0.5) == True:
                result.append(str(cytoband_chr)+cytoband_info)
        except:
            print("ERROR")
            pdb.set_trace()

    return result

def fetchReadDepth(bam_cram_file, t_chr, t_start, t_stop, MAQ):
    name, ext = os.path.splitext(bam_cram_file)
    if ext == '.cram':
        f = pysam.AlignmentFile(bam_cram_file, 'rc')
    elif exit == '.bam':
        f = pysam.AlignmentFile(bam_cram_file, 'rb')

    if not f.has_index():
        print('No index found for ', bam_cram_file)
        sys.exit(0)
    try:
        rd_list = []
        num = 0
## Method1 fetch read count. Not match CN-Learn results
#        iter = f.fetch(t_chr, t_start, t_stop)
#        for i in iter:
#            if i.maq >= MAQ:
#                num += 1
## Method2 fetch read depth. Not match CN-Learn results
#        for i in f.pileup(t_chr, t_start, t_stop,until_eof=True):
#            if i.pos >= t_start and i.pos <= t_stop:
#                num += 1
#                print ("%d: Coverage at base %s = %s" % (num, i.pos, i.n))
#                pdb.set_trace()
#                rd_list.append(i.n)
#
## Method3 pysam samtools
        for pileup_pos in range(t_start+1,t_stop+1):
            num += 1
            rd = pysam.view(bam_cram_file,t_chr+":"+str(pileup_pos)+"-"+str(pileup_pos))
            rd = rd.strip()
            rd = len(rd.split("\n"))
            rd_list.append(rd)
    except:
        print("[ERROR] Could not retrieve mappings for region %s:%d-%d.\
                Check that contigs are named correctly and the bam file \
                is properly indexed" % (t_chr,t_start,t_stop))
        sys.exit(0)
    mean_rd = np.average(rd_list)
    return mean_rd

def fetch_consolidated_calls(sampleID, cnv_chr, cnv_start, cnv_stop,cnv_type, consolidated_cnv_list):
    # function: fetch cn_learn consolidated calls
    # 2020.6.15
    overlap_w_canoes,overlap_w_clamms,overlap_w_xhmm = 0,0,0
    for cons_cnv_reader in consolidated_cnv_list:
        cons_cnv_chr = cons_cnv_reader[0]
        cons_cnv_start = int(cons_cnv_reader[1])
        cons_cnv_stop = int(cons_cnv_reader[2])
        cons_cnv_type = cons_cnv_reader[3]
        cons_cnv_caller = cons_cnv_reader[5]

        overlap_region = reciprocalOverlap_returnOverlap(cnv_chr, cnv_start, cnv_stop, cnv_type, cons_cnv_chr, cons_cnv_start, 
                cons_cnv_stop, cons_cnv_type, 1)
        if overlap_region != 0:
            # Note: the ratio is divied by the length of CNV prediced by CN-Learn, not each tools' prediction
            # so using "overlap_ratio = round(float(overlap_region/(cons_cnv_stop-cons_cnv_start)),2)" is problematic
            overlap_ratio = round(float(overlap_region/(cnv_stop-cnv_start)),2)
            if overlap_ratio == 1.0:
                overlap_ratio = 1

            if cons_cnv_caller == "XHMM":
                overlap_w_xhmm = overlap_ratio 
            elif cons_cnv_caller == "CANOES":
                overlap_w_canoes = overlap_ratio 
            elif cons_cnv_caller == "CLAMMS":
                overlap_w_clamms = overlap_ratio 
            else:
                pdb.set_trace()
    return [overlap_w_canoes,overlap_w_clamms,overlap_w_xhmm]

def fetch_targets_num(cnv_chr, cnv_start, cnv_stop, target_probe_list):
    num = -1 #list start from 0
    cnv_target_num_list = []
    target_start_num, target_stop_num = None, None
    for target_reader in target_probe_list:
        num += 1
        target_chr = target_reader[0]
        target_start = int(target_reader[1])
        target_stop = int(target_reader[2])
    #     if cnv_chr == target_chr:
    #         # if target_start in range(cnv_start, cnv_stop+1) and target_stop in range(cnv_start, cnv_stop+1):
    #         #     cnv_target_num_list.append(num) #for CNVs span many targets
    #         # elif cnv_start in range(target_start, target_stop+1) and cnv_stop in range(target_start, target_stop+1):
    #         #     cnv_target_num_list.append(num) #for one target CNVs

    #         # Three scenarios: 1) cnv span one target; 2) cnv span two targets; and 3)cnv span may targets.
    #         # Pls note that CLAMMS divides exome capture regions that are ≥ 1000 bp long into equally-sized 500-1000 bp windows. 
    #         # therefore the cnv and target probe boundaries may be not same.
    #         if cnv_start in range(target_start, target_stop+1) or cnv_stop in range(target_start, target_stop+1):
    #             pdb.set_trace()
    #             cnv_target_num_list.append(num) #for CNVs span many targets
    #         else:
    #             pass
    #     else:
    #         if cnv_target_num_list != []:
    #             break
        if reciprocalOverlap(cnv_chr, cnv_start, cnv_stop, None, target_chr, target_start, target_stop, None, 1):
            cnv_target_num_list.append(num)

    if cnv_target_num_list == []:
        pdb.set_trace()
    return cnv_target_num_list

def fetch_targets_num_w_flanking(cnv_target_num_list, adjust_target_num):
    left_flank_target_num_list,right_flank_target_num_list=[],[]
    try:
        left_flank_target_num_list = list(range(cnv_target_num_list[0]-adjust_target_num, cnv_target_num_list[0]))
        right_flank_target_num_list = list(range(cnv_target_num_list[-1]+1, cnv_target_num_list[-1]+adjust_target_num+1))
        #TODO: what should return if the CNV start/stop target is the first/last target in a chromosome?
    except:
        pdb.set_trace()
    return [left_flank_target_num_list,right_flank_target_num_list]

def fetchReadDepthFromTabularFile(tabular_file, t_chr, t_start, t_stop):
    if not os.path.exists(tabular_file):
        print('No tabular file: %s'%tabular_file)
        sys.exit(0)
    if not os.path.exists(tabular_file+'.tbi'):
        pysam.tabix_index(tabular_file,seq_col=0, start_col=1, end_col=2)
    f = pysam.TabixFile(tabular_file)
    rd_list = []
    num = 0
    for row in f.fetch(t_chr, t_start, t_stop, parser=pysam.asTuple()):
        row_start = int(row[1])
        row_stop = int(row[2])
        for iter in range(max(row_start,t_start), min(row_stop,t_stop)):
            rd_list.append(int(row[3]))
    return np.mean(rd_list)    

def fetch_targets_rd(bp_coverage_file, target_num_list, target_probe_list):
    rd_list = []
    for item in target_num_list:
        try:
            target_reader = target_probe_list[item]
            target_chr = target_reader[0]
            target_start = int(target_reader[1])
            target_stop = int(target_reader[2])
            target_rd = fetchReadDepthFromTabularFile(bp_coverage_file, target_chr, target_start, target_stop)
            rd_list.append(target_rd)
        except:
            pass
        
    rd = np.mean(rd_list)
    return rd

def fetch_rd_ratio(bp_coverage_file, cnv_chr, cnv_start, cnv_stop, num_targets, target_probe_list):
    if os.path.isfile(bp_coverage_file):
        half_num_target = max(int(num_targets/2),1) #TODO: need to check the results with CN-learn
        cnv_target_num_list = fetch_targets_num(cnv_chr, cnv_start, cnv_stop, target_probe_list)
        left_flank_target_num_list,right_flank_target_num_list = fetch_targets_num_w_flanking(cnv_target_num_list, half_num_target)
        cnv_target_rd = fetch_targets_rd(bp_coverage_file, cnv_target_num_list, target_probe_list)
        left_flank_rd = fetch_targets_rd(bp_coverage_file, left_flank_target_num_list, target_probe_list)
        right_flank_rd = fetch_targets_rd(bp_coverage_file, right_flank_target_num_list, target_probe_list)
        rd_ratio = round(cnv_target_rd/((left_flank_rd+right_flank_rd)/2),2)
    else:
        rd_ratio = None

    return rd_ratio

def fetch_rc_ratio(rc_ratio_file, cnv_chr, cnv_start, cnv_stop):
    rc_ratio_list = []
    if not os.path.exists(rc_ratio_file):
        print('No tabular file: %s'%rc_ratio_file)
        return ["No tabular file"]
    if not os.path.exists(rc_ratio_file+'.tbi'):
        pysam.tabix_index(rc_ratio_file,seq_col=0, start_col=1, end_col=2)
    f = pysam.TabixFile(rc_ratio_file)
    rc_ratio_list = []
    for row in f.fetch(cnv_chr, int(cnv_start), int(cnv_stop), parser=pysam.asTuple()):
        rc_ratio_list.append(float(row[3]))
    return rc_ratio_list


# function
def sampleHasVarints(sample_ID, vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file,'r'))
    samples_list = vcf_reader.samples
    if sample_ID in samples_list:
        return True
    else:
        return False
       
def event_stats(xhmm_vcfdict_cnv, Q):
    cnv_sq = float(xhmm_vcfdict_cnv[5])
    cnv_ndq = float(xhmm_vcfdict_cnv[8])

    if (cnv_sq >= Q) and (cnv_ndq < Q):
        status = "HAS_EVENT"
    elif not (cnv_sq < Q) and (cnv_ndq >= Q):
        status = "DOES_NOT_HAVE_EVENT"
    else:
        status = "MISSING"
    return status

def reciprocalOverlap(chrom_data1, start_data1, stop_data1, type_data1, chrom_data2, start_data2, stop_data2, type_data2, thresthold):
    if chrom_data1 == chrom_data2:
        start_data1 = int(start_data1)
        stop_data1 = int(stop_data1)
        start_data2 = int(start_data2)
        stop_data2 = int(stop_data2)
        
        if type_data1 == type_data2: 
            if start_data1 in range(start_data2, stop_data2+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data1 + 1
                if thresthold == 1:  #1bp overlap
                    if overlap > 0:
                        return True
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        #print overlap/(stop_data1-start_data1+1),overlap/(stop_data2-start_data2+1)
                        return True
            elif start_data2 in range(start_data1, stop_data1+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data2 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return True
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            elif start_data1 in range(start_data2, stop_data2+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data1 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return True
                else:    
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            elif start_data2 in range(start_data1, stop_data1+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data2 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return True
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            else:
                return False
        else:
            return False
    else:
        return False

def reciprocalOverlap_returnOverlap(chrom_data1, start_data1, stop_data1, type_data1, chrom_data2, start_data2, stop_data2, type_data2, thresthold):
    if chrom_data1 == chrom_data2:
        start_data1 = int(start_data1)
        stop_data1 = int(stop_data1)
        start_data2 = int(start_data2)
        stop_data2 = int(stop_data2)
        if type_data1 == type_data2: 
            if start_data1 in range(start_data2, stop_data2+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data1 
                if thresthold == 1:  #1bp overlap
                    if overlap > 0:
                        return overlap
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return overlap
            elif start_data2 in range(start_data1, stop_data1+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data2
                if thresthold == 1:
                    if overlap > 0:
                        return overlap
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return True
            elif start_data1 in range(start_data2, stop_data2+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data1 
                if thresthold == 1:
                    if overlap > 0:
                        return overlap
                else:    
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return overlap
            elif start_data2 in range(start_data1, stop_data1+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data2
                if thresthold == 1:
                    if overlap > 0:
                        return overlap
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) and (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return overlap
            else:
                return 0
        else:
            return 0
    else:
        return 0
        
def regionOverlap_occupyFirst(chrom_data1, start_data1, stop_data1, chrom_data2, start_data2, stop_data2, thresthold):
    if chrom_data1 == chrom_data2:
        start_data1 = int(start_data1)
        stop_data1 = int(stop_data1)
        start_data2 = int(start_data2)
        stop_data2 = int(stop_data2)

        if start_data1 in range(start_data2, stop_data2+1) and stop_data2 in range(start_data1, stop_data1+1):
            overlap = stop_data2 - start_data1 + 1
            if thresthold == 1:  #1bp overlap
                if overlap > 0:
                    return True
            else:
                if (overlap/(stop_data1-start_data1+1) >= thresthold): # or (overlap/(stop_data2-start_data2+1)>=thresthold)
                    return True
        elif start_data2 in range(start_data1, stop_data1+1) and stop_data1 in range(start_data2, stop_data2+1):
            overlap = stop_data1 - start_data2 + 1
            if thresthold == 1:
                if overlap > 0:
                    return True
            else:
                if (overlap/(stop_data1-start_data1+1) >= thresthold): # or (overlap/(stop_data2-start_data2+1)>=thresthold)
                    return True
        elif start_data1 in range(start_data2, stop_data2+1) and stop_data1 in range(start_data2, stop_data2+1):
            overlap = stop_data1 - start_data1 + 1
            if thresthold == 1:
                if overlap > 0:
                    return True
            else:    
                if (overlap/(stop_data1-start_data1+1) >= thresthold): # or (overlap/(stop_data2-start_data2+1)>=thresthold)
                    return True
        elif start_data2 in range(start_data1, stop_data1+1) and stop_data2 in range(start_data1, stop_data1+1):
            overlap = stop_data2 - start_data2 + 1
            if thresthold == 1:
                if overlap > 0:
                    return True
            else:
                if (overlap/(stop_data1-start_data1+1) >= thresthold): # or (overlap/(stop_data2-start_data2+1)>=thresthold)
                    return True
        else:
            return False
    else:
        return False

def overlap_condition(chrom_data1, start_data1, stop_data1, chrom_data2, start_data2, stop_data2, thresthold):
    if chrom_data1 == chrom_data2:
        try:
            if start_data1 in range(start_data2, stop_data2+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data1 + 1
                if thresthold == 1:  #1bp overlap
                    if overlap > 0:
                        return "overlap_3UTR"
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) or (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return "overlap_3UTR"

            elif start_data2 in range(start_data1, stop_data1+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data2 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return "overlap_5UTR"
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) or (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return "overlap_5UTR"

            elif start_data1 in range(start_data2, stop_data2+1) and stop_data1 in range(start_data2, stop_data2+1):
                overlap = stop_data1 - start_data1 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return "overlap_middle"
                else:    
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) or (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return "overlap_middle"

            elif start_data2 in range(start_data1, stop_data1+1) and stop_data2 in range(start_data1, stop_data1+1):
                overlap = stop_data2 - start_data2 + 1
                if thresthold == 1:
                    if overlap > 0:
                        return "overlap_complete"
                else:
                    if (overlap/(stop_data1-start_data1+1) >= thresthold) or (overlap/(stop_data2-start_data2+1)>=thresthold):
                        return "overlap_complete"
            
            elif (start_data1 == start_data2) and (stop_data1 == stop_data2):
                return "overlap_identical_complete"
                
            else:
                return "-"

        except Exception as e:
            pdb.set_trace()
            raise e
        
    else:
        return False

def cnv_length(cnv_row):
    start = cnv_row[2]
    stop = cnv_row[3]
    cnv_length = int(stop-start+1)
    return cnv_length

def find_denovo_largeInherited_CNV(family_ID, child_cnv_list, father_cnv_list, mother_cnv_list, o_threshold, l_threshold, Q):
    denovo_cnv_results = []
    large_inherited_cnv_results = []
    # kilobases
    l_threshold = l_threshold*1000
    print("Overlap Threshold:",o_threshold, "\tLength Thresthold:", l_threshold)
    for child_cnv in child_cnv_list:
        denovo_flag = True
        for mother_cnv in mother_cnv_list:
            if reciprocal_length(child_cnv,mother_cnv,o_threshold) == True:
                denovo_flag = False
                break
        if denovo_flag == True:
            for father_cnv in father_cnv_list:
                if reciprocal_length(child_cnv,father_cnv,o_threshold) == True:
                    denovo_flag = False
                    break

        family_list = []
        family_list.append(family_ID)
        family_list.extend(child_cnv)
        if denovo_flag == True:
            denovo_cnv_results.append(family_list)
        elif denovo_flag == False and cnv_length(child_cnv) >= l_threshold:
            large_inherited_cnv_results.append(family_list)  
    return [denovo_cnv_results, large_inherited_cnv_results]

def isParent(ped_dict, sample_id):
    if sample_id in ped_dict:
        pedigree_info = ped_dict[sample_id]
        if pedigree_info[2] == '0' and pedigree_info[3] == '0' and pedigree_info[6] == 'Parent':
            return True
        else:
            return False
    else:
        pass

def isOffspring(ped_dict, sample_id):
    if sample_id in ped_dict:
        pedigree_info = ped_dict[sample_id]
        if pedigree_info[6] == 'Offspring':
            return True
        else:
            return False
    else:
        return None

def isAffected(ped_dict, sample_id):
    if sample_id in ped_dict:
        pedigree_info = ped_dict[sample_id]
        phenotype = pedigree_info[5]
        if phenotype == '1':
            return "unaffected"
        elif phenotype == '2':
            return "affected"
        elif phenotype == '0':
            return "missing"
        elif phenotype == '-9':
            return "missing"
        else:
            return phenotype
    else:
        return None

def buildXHMMrdDictIndex(XHMMrd_dict):
    coordinate_list = XHMMrd_dict[b'Matrix']
    pattern = re.compile(rb'(\S+):(\d+)-(\d+)')
    num = 0
    target_list=[]
    for target_reader in coordinate_list:
        num += 1
        target = pattern.match(target_reader)
        exon_chr = target.group(1)
        exon_start = int(target.group(2))
        exon_stop = int(target.group(3))
