import vcf
import pdb

def min(a, b):
        return a if a <= b else b

def max(a, b):
    return a if a >= b else b

def chrStr2Int(chr):
    chr = chr.replace('chr','')
    if chr == 'X' or chr == 'x':
        return 23
    elif chr == 'Y' or chr == 'y':
        return 24
    else:
        return int(chr)

class VCFReader:
    def __init__(self, vcf_file):
        try:
            vfile = open(vcf_file, 'r')
            self.vcf_reader = vcf.Reader(vfile)
        except IOError:
            print 'This VCF file does not exist: ' + vcf_file

    def chrInt2Str(self, chromosome_int):
        try:
            if int(chromosome_int) == 23:
                return 'X'
            elif int(chromosome_int) == 24:
                return 'Y' 
            else:
                return str(chromosome_int)
        except:
            return str(chromosome_int)
        

    # Given a region, check the number of sites which are heterozygous calls
    def getSNPInfo(self, sample_id, intervals_list):
        snp_info = list()
        if self.vcf_reader:
            for in_list in intervals_list:
                for interval in in_list:
                    # The start and end coordinates are in the zero-based, half-open coordinate system
                    try:
                        records = self.vcf_reader.fetch(self.chrInt2Str(interval._chr), interval._bp1, interval._bp2 + 1) 
                        if records:
                            het_num = 0
                            for record in records:
                                genotype = record.genotype(sample_id)
                                if genotype and genotype.is_het:
                                    het_num += 1    
                            snp_info.append(het_num)
                        else:
                            snp_info.append(0)
                    except Exception, e:
                        print 'WARRNING:',e                    
                        
        return snp_info
    
    def loadTagSNP(self, filename):
        tagsnp_f=open(filename,'r')
        cnp_dict = {}

        next(tagsnp_f) # Skip first row
        for row in tagsnp_f:
            row=row.strip()
            element = row.split('\t')
            cnp_name = element[0]
            chr = element[1].replace("chr","")
            cnp_start = int(element[2])
            cnp_stop = int(element[3])
            snp_pos = int(element[4])
            R2 = float(element[5])
            population = element[6]

            if cnp_dict.has_key(population):
                cnp_dict[population].append([cnp_name,chr,cnp_start,cnp_stop,snp_pos,R2,population])
            else:
                cnp_dict[population] = [[cnp_name,chr,cnp_start,cnp_stop,snp_pos,R2,population]]

        return cnp_dict

    def findTagSNPForSample(self, sample_p, sample_id, cnp_dict):
        #Renjie modified
        cnp_list = []
        
        if self.vcf_reader:
            for cnp in cnp_dict[sample_p]:
                records = self.vcf_reader.fetch(cnp[1], cnp[4]-1, cnp[4]) 
                if records:
                    for record in records:
                        #Renjie added try ... except
                        try:
                            genotype = record.genotype(sample_id)['GT']
                            if not genotype == '0|0' and not genotype == '0/0':
                                cnp_list.append(cnp) 
                        except Exception, e:
                            print 'WARRNING:',e   

        return cnp_list 


    def findExonWithTagSNP(self, cnp_list, targets_list, overlap_threshold):
        overlap_exon = []
        overlap_exon_unique = []
        
        tagsnp_result = {} 

        chr_list = []
        chr_targets = {}
        for targets in targets_list:
            tar_chr = chrStr2Int(targets[0]._chr)
            chr_list.append(tar_chr)
            chr_targets[tar_chr] = targets
            tagsnp_result[tar_chr] = [0.0] * len(targets)

            resuts = [[0.0]] * len(targets)
        for cnp in cnp_list:
            cnp_name = cnp[0]
            cnp_chr = chrStr2Int(cnp[1])
            cnp_start = cnp[2]
            cnp_stop = cnp[3]
            snp_pos = cnp[4]
            R2 = float(cnp[5])
            population = cnp[6]

            try:
                targets = chr_targets[cnp_chr]
                results = tagsnp_result[cnp_chr]
            except KeyError:
                print 'Warning: No chromosome %s in exome_region_dict' %cnp_chr

            for i, exon in enumerate(targets):
                exon_chr = exon._chr
                exon_start = exon._bp1 + 1 #for bed file
                exon_stop = exon._bp2 + 1 #for bed file
                exon_size = exon_stop - exon_start + 1
                
                if exon_start > cnp_stop:
                    break
                if exon_stop >= cnp_start:
                    o_start = max(exon_start, cnp_start)
                    o_stop = min(exon_stop, cnp_stop)

                    overlap_ratio = float(o_stop - o_start + 1) / float(exon_size)

                    if overlap_ratio >= overlap_threshold:
                        if R2 >= 0.8 and R2 > results[i]:
                            results[i] = R2 #TODO if one exon is overlapped twice, how to merge several R2 values ?
                            temp = cnp[:]
                            e_list = [exon._chr, exon._bp1, exon._bp2]
                            temp.extend(e_list)
                            overlap_exon.append(temp)
                            #print temp
                            #store the unique data
                            if e_list not in overlap_exon_unique:
                                overlap_exon_unique.append(e_list)
        tagsnp_info = []
        for chr in chr_list:
            tagsnp_info.append(tagsnp_result[chr])

        return tagsnp_info
