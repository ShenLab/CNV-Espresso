from __future__ import division
import os
import argparse
import glob
import function as df
from DataManager import *
from hmm.Model import *
from hmm.ModelParams import *
import operator 
import numpy as np
from VCFReader import *
from ParameterEstimation import *
import fileinput
import pdb
import sys
import time

def bamlist2RPKM(args):
    MAQ = args.maq 
    print "MAQ threshold:",MAQ

    try:
        import pysam
    except:
        print 'Cannot load pysam module!'
        sys.exit(0)
    try:
        # read target
	target_fn = str(args.target)
	targets = jf.loadTargets(target_fn)
	num_target = len(targets)
    except IOError:
        print 'Cannot read target file: ', target_fn
        sys.exit(0)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    try:
        bamlist_f = open(args.input)
    except IOError:
        sys.exit(0)

    for line in bamlist_f.readlines():
        line = line.strip('\n')
        temp = line.split('\t')
        sample_name = temp[0] + '.' + temp[2]
        bam_file = line.split('\t')[1]
        
        f = pysam.AlignmentFile(bam_file, 'rb')

        if not f.has_index():
            print 'No index found for ', bam_file
            sys.exit(0)
    
        readcount = np.zeros(num_target)
        exon_bp = np.zeros(num_target)
        targetIDs = np.zeros(num_target)

        # detect contig naming scheme here # TODO, add an optional "contigs.txt" file or automatically handle contig naming
        bam_contigs = f.references
        targets_contigs = [str(t) for t in set(map(operator.itemgetter("chr"), targets))]
                    
        targets2contigmap = {}
        for targets_contig in targets_contigs:
            if targets_contig in bam_contigs:
                targets2contigmap[targets_contig] = targets_contig
            elif jf.chrInt2Str(targets_contig) in bam_contigs:
                targets2contigmap[targets_contig] = jf.chrInt2Str(targets_contig)
            elif jf.chrInt2Str(targets_contig).replace("chr","") in bam_contigs:
                targets2contigmap[targets_contig] = jf.chrInt2Str(targets_contig).replace("chr","")
            else:
                print "[ERROR] Could not find contig '%s' from %s in bam file! \n[ERROR] Perhaps the contig names for the targets are incompatible with the bam file ('chr1' vs. '1'), or unsupported contig naming is used?" % (targets_contig, target_fn)
                sys.exit(0)

        print 'Calculating RC and RPKM values...'
        
        total_reads = 0
        counter = 0
        for t in targets:
            t_chr = targets2contigmap[str(t['chr'])]

            t_start = t['start']
            t_stop = t['stop']
            
            try:
                iter = f.fetch(t_chr, t_start, t_stop)
            except:
                print "[ERROR] Could not retrieve mappings for region %s:%d-%d. Check that contigs are named correctly and the bam file is properly indexed" % (t_chr,t_start,t_stop)
                sys.exit(0)

            for i in iter:
                if i.pos+1 >= t_start and i.mapq >= MAQ:
                    readcount[counter] += 1
                    total_reads += 1

            exon_bp[counter] = t_stop - t_start + 1
            targetIDs[counter] = counter + 1
            counter += 1

        print 'Found %d reads in the target regions of bam file with MAQ >= %d: ' %(total_reads, MAQ), bam_file
        # calculate RPKM values for all targets 
        rpkm = (10**9*(readcount)/(exon_bp))/(total_reads)
        rpkm_f = open(args.output+'/'+sample_name+'.rc.rpkm', 'w')
        
        rpkm_f.write('chr\tstart\tstop\tRC\tRPKM\n')
        for i in range(len(rpkm)):
            rpkm_f.write(targets[i]['chr'] + '\t' + str(targets[i]['start']) + '\t' + str(targets[i]['stop']) + '\t' + str(readcount[i]) + '\t' + str(rpkm[i]) + '\n')
        rpkm_f.close()

    bamlist_f.close()
    
def RPKM2Matrix(args):
    #Renjie modified
    rpkm_dir = str(args.rpkm_dir)
    rpkm_files = glob.glob(rpkm_dir + "/*")
    if len(rpkm_files) == 0:
        print 'Cannot find any rpkm files'
        sys.exit(0)
    if not os.path.exists(args.output):
        os.mkdir(args.output)
        print 'Output dir created: ',args.output

    try:
        # read target
        target_fn = str(args.target)
        targets = jf.loadTargetsStr(target_fn)
        num_target = len(targets)
    except IOError:
        print 'Cannot read target file: ', target_fn
        sys.exit(0)
   
    samples = {}
    for f in rpkm_files:
        s = '.'.join(f.split('/')[-1].split('.')[0:-1])
        samples[s] = f

    RPKM_matrix = np.zeros([num_target, len(samples)], dtype=np.float)
    RC_matrix = np.zeros([num_target, len(samples)], dtype=np.float)
    for i,s in enumerate(samples.keys()):
        rc = np.loadtxt(samples[s], dtype=np.float, delimiter="\t", skiprows=1, usecols=[3])
        rpkm = np.loadtxt(samples[s], dtype=np.float, delimiter="\t", skiprows=1, usecols=[4])
        RC_matrix[:,i] = rc
        RPKM_matrix[:,i] = rpkm
        print "Successfully read RC and RPKM for " + s
    
    output_rpkm = '/RPKM_matrix.rpkm.raw'
    output_rc = '/RC_matrix.rc.raw'

    if args.output:
        output_rpkm = args.output+'/RPKM_matrix.raw'
        output_rc = args.output+'/RC_matrix.raw'
    output_rc_f = open(output_rc, 'w')
    output_rpkm_f = open(output_rpkm, 'w')

    output_rc_f.write('Targets\t' + '\t'.join(samples.keys()) + '\n')
    for i in range(len(RC_matrix)):
        output_rc_f.write(targets[i] + '\t' + '\t'.join(str(r) for r in RC_matrix[i]) + '\n')
    output_rc_f.close()

    output_rpkm_f.write('Targets\t' + '\t'.join(samples.keys()) + '\n')
    for i in range(len(RPKM_matrix)):
        output_rpkm_f.write(targets[i] + '\t' + '\t'.join(str(r) for r in RPKM_matrix[i]) + '\n')
    output_rpkm_f.close()
    
def filter_rpkm(args):
    #Renjie revised
    rpkm_matrix = str(args.rpkm_matrix)
    raw_dir = os.path.dirname(rpkm_matrix)
    if raw_dir != '':
        raw_dir = raw_dir + '/'

    output = rpkm_matrix
    if args.output:
        output = str(args.output)
    print 'Loading targets...'
    temp = jf.loadTargetsFromFirstCol(rpkm_matrix)
    targets = temp['targets']
    targets_str = temp['targets_str']

    print 'Loading parameters from ' + args.filter_params
    f = open(args.filter_params)
    params = f.readline().split('\t')
    min_gc = int(params[0])
    max_gc = int(params[1])
    min_map = int(params[2])
    max_exon = int(params[3])
    min_rpkm = float(params[4])

    #Fileter targets which median of RPKM < min_rpkm
    print "Calculating GC content..."
    try:
        GC_percentage = jf.loadNormValues(raw_dir + 'GC_percentage')
    except:
        GC_percentage = jf.calGCPercentage(targets, args.ref_file)
        jf.saveNormValues(raw_dir + 'GC_percentage', targets_str, GC_percentage, 'GC_content')
    
    print 'Calculating mapping ability...'
    try:
        map_ability = jf.loadNormValues(raw_dir + 'mapping_ability')
    except:
        map_ability = jf.calMapAbility(targets, args.map_file)
        jf.saveNormValues(raw_dir + 'mapping_ability', targets_str, map_ability, 'Mapping_ability')
    
    print 'Calculating exon length...'
    try:
        exon_length = jf.loadNormValues(raw_dir + 'exon_length')
    except:
        exon_length = jf.calExonLength(targets)
        jf.saveNormValues(raw_dir + 'exon_length', targets_str, exon_length, 'Exon_length')

    print 'Filtering targets by GC content, mapping ability and exon length'
    matrix = open(rpkm_matrix)
    n_lines = []
    e_lines = []
    lines = matrix.readlines()
    t = lines[0].index('\t')
    title = lines[0][0:t] + '\t' + '\t'.join(['GC_Content', 'Mapping_Ability', 'Exon_Length']) + lines[0][t:]
    for n in range(len(lines)-1):
        line = lines[n+1]
        t = line.index('\t')
        gc = GC_percentage[n]
        m = map_ability[n]
        el = exon_length[n]
        rpkm_line = line.strip().split('\t')[1:]
        rpkm_line_median = np.median([float(s) for s in rpkm_line])
        l = line[0:t] + '\t' + '\t'.join([str(gc), str(m), str(el)]) + line[t:]
        if gc < min_gc or gc > max_gc:
            l = l.strip() + '\t' + 'RemovedByGC' +'\n'
            e_lines.append(l)
        elif m < min_map:
            l = l.strip() + '\t' + 'RemovedByMAP'+'\n'
            e_lines.append(l)
        elif el > max_exon:
            l = l.strip() + '\t' + 'RemovedByExon_length'+'\n'
            e_lines.append(l)
        elif rpkm_line_median < min_rpkm:
            l = l.strip() + '\t' + 'RemovedByMedian_rpkm'+'\n'
            e_lines.append(l)
        else:
            n_lines.append(l)

    n_file = open(output + '.filtered', 'w')
    n_file.write(title)
    n_file.writelines(n_lines)
    n_file.close()

    e_file = open(output + '.outliers', 'w')
    e_file.write(title)
    e_file.writelines(e_lines)
    e_file.close()
                
def normalize(args): 
    rpkm_matrix = str(args.rpkm_matrix)
    raw_dir = os.path.dirname(rpkm_matrix)
    if raw_dir != '':
        raw_dir = raw_dir + '/'

    output = rpkm_matrix + '.normalized'
    if args.output:
        output = raw_dir + str(args.output) + '.normalized'

    print 'Loading matrix...'
    result = jf.loadRPKMMatrix(rpkm_matrix)
    rpkm = result['rpkm']
    samples = result['samples']
    annotation = result['annotation']
    targets = result['targets']
    targets_str = result['targets_str']

    GC_percentage = annotation[:, 0]
    map_ability = annotation[:, 1]
    exon_length = annotation[:, 2]

    GC_index = {}
    for ind in range(len(GC_percentage)):
        gc = GC_percentage[ind]

        if GC_index.has_key(gc):
            GC_index[gc].append(ind)
        else:
            GC_index[gc] = [ind]

    print 'Normalizing by GC percentage...'
    corrected_rpkm = np.zeros([len(rpkm), len(rpkm[0])], dtype=np.float)
    for i in range(len(samples)):
        print 'Normalizing RPKM by GC content for sample: ' + samples[i]
        overall_median = np.median(rpkm[:, i])
        for gc in GC_index.keys():
            t_ind = GC_index[gc]
            t_median = np.median(rpkm[t_ind, i])
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, GC: %d' %(samples[i], gc)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = rpkm[t_ind, i] * overall_median / t_median

    print 'Saving GC normalized matrix..'
    file_GC_normalized = open(output + '.GC', 'w')
    np.savetxt(file_GC_normalized, corrected_rpkm, delimiter='\t', header='\t'.join(samples),comments='')
    file_GC_normalized.close()


    map_index = {}
    for ind in range(len(map_ability)):
        _map = map_ability[ind]

        if map_index.has_key(_map):
            map_index[_map].append(ind)
        else:
            map_index[_map] = [ind]

    print 'Normalizing by Mapping ability...'
    for i in range(len(samples)):
        print 'Normalizing RPKM by mapping ability for sample %s' %samples[i]
        overall_median = np.median(corrected_rpkm[:, i])
        for _map in map_index.keys():
            t_ind = map_index[_map]
            t_median = np.median(corrected_rpkm[t_ind, i])
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, Mapping ability: %d' %(samples[i], _map)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = corrected_rpkm[t_ind, i] * overall_median / t_median

    print 'Saving Mappability normalized matrix..'
    file_MAP_normalized = open(output + '.MAP', 'w')
    np.savetxt(file_MAP_normalized, corrected_rpkm, delimiter='\t', header='\t'.join(samples),comments='')
    file_MAP_normalized.close()

    length_index = {}
    for ind in range(len(exon_length)):
        _length = exon_length[ind]
        if length_index.has_key(_length):
            length_index[_length].append(ind)
        else:
            length_index[_length] = [ind]

    print 'Normalizing by exon length...'
    for i in range(len(samples)):
        print 'Normalizing RPKM for by exon length for sample %s' %samples[i]
        overall_median = np.median(corrected_rpkm[:, i])
        for _length in length_index.keys():
            t_ind = length_index[_length]
            t_median = np.median(corrected_rpkm[t_ind, i])
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, Exome length: %d' %(samples[i], _length)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = corrected_rpkm[t_ind, i] * overall_median / t_median
    
    print 'Saving exon_length normalized matrix..'
    file_exon_length_normalized = open(output + '.exon_length', 'w')
    np.savetxt(file_exon_length_normalized, corrected_rpkm, delimiter='\t', header='\t'.join(samples),comments='')
    file_exon_length_normalized.close()

    print 'Saving matrix..'
    jf.saveRPKMMatrix(output, samples, targets_str, np.transpose(corrected_rpkm))
    

def normalize_gc(args): 
    rpkm_matrix = str(args.rpkm_matrix)
    raw_dir = os.path.dirname(rpkm_matrix)
    if raw_dir != '':
        raw_dir = raw_dir + '/'

    output = rpkm_matrix + '.normalized'
    if args.output:
        output = raw_dir + str(args.output) + '.normalized'

    print 'Loading matrix...'
    result = jf.loadRPKMMatrix(rpkm_matrix)
    rpkm = result['rpkm']
    samples = result['samples']
    annotation = result['annotation']
    targets = result['targets']
    targets_str = result['targets_str']

    GC_percentage = annotation[:, 0]
    map_ability = annotation[:, 1]
    exon_length = annotation[:, 2]

    GC_index = {}
    for ind in range(len(GC_percentage)):
        gc = GC_percentage[ind]

        if GC_index.has_key(gc):
            GC_index[gc].append(ind)
        else:
            GC_index[gc] = [ind]

    print 'Normalizing by GC percentage...'
    corrected_rpkm = np.zeros([len(rpkm), len(rpkm[0])], dtype=np.float)
    for i in range(len(samples)):
        print 'Normalizing RPKM by GC content for sample: ' + samples[i]
        overall_median = np.median(rpkm[:, i])
        for gc in GC_index.keys():
            t_ind = GC_index[gc]
            t_median = np.median(rpkm[t_ind, i])
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, GC: %d' %(samples[i], gc)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = rpkm[t_ind, i] * overall_median / t_median

    print 'Saving GC normalized matrix..'
    file_GC_normalized = open(output + '.GC', 'w')
    np.savetxt(file_GC_normalized, corrected_rpkm, delimiter='\t', header='\t'.join(samples),comments='')
    file_GC_normalized.close()


    map_index = {}
    for ind in range(len(map_ability)):
        _map = map_ability[ind]

        if map_index.has_key(_map):
            map_index[_map].append(ind)
        else:
            map_index[_map] = [ind]

    print 'Normalizing by Mapping ability...'
    for i in range(len(samples)):
        print 'Normalizing RPKM by mapping ability for sample %s' %samples[i]
        overall_median = np.median(corrected_rpkm[:, i])
        for _map in map_index.keys():
            t_ind = map_index[_map]
            t_median = np.median(corrected_rpkm[t_ind, i])
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, Mapping ability: %d' %(samples[i], _map)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = corrected_rpkm[t_ind, i] * overall_median / t_median

    print 'Saving Mappability normalized matrix..'
    file_MAP_normalized = open(output + '.MAP', 'w')
    np.savetxt(file_MAP_normalized, corrected_rpkm, delimiter='\t', header='\t'.join(samples),comments='')
    file_MAP_normalized.close()

    length_index = {}
    for ind in range(len(exon_length)):
        _length = exon_length[ind]
        if length_index.has_key(_length):
            length_index[_length].append(ind)
        else:
            length_index[_length] = [ind]

    print 'Normalizing by exon length...'
    for i in range(len(samples)):
        print 'Normalizing RPKM for by exon length for sample %s' %samples[i]
        overall_median = np.median(corrected_rpkm[:, i])
        for _length in length_index.keys():
            t_ind = length_index[_length]
            t_median = np.median(corrected_rpkm[t_ind, i])
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, Exome length: %d' %(samples[i], _length)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = corrected_rpkm[t_ind, i] * overall_median / t_median
    
    print 'Saving exon_length normalized matrix..'
    file_exon_length_normalized = open(output + '.exon_length', 'w')
    np.savetxt(file_exon_length_normalized, corrected_rpkm, delimiter='\t', header='\t'.join(samples),comments='')
    file_exon_length_normalized.close()

    print 'Saving matrix..'
    jf.saveRPKMMatrix(output, samples, targets_str, np.transpose(corrected_rpkm))
    
def svd(args):
    filename = args.rpkm_matrix 
    f_dir = os.path.dirname(filename)
    if f_dir != '':
        f_dir = f_dir + '/'

    output = filename
    if args.output:
        output = f_dir + str(args.output)
   
    # count the columns number of the data file
    f = open(filename)
    temp = f.readline().strip().split('\t')
    colsnum = len(temp)
    
    # skip 1st row and 4 columns
    print 'Loading file...'
    data = np.loadtxt(filename, dtype=np.float, delimiter='\t', skiprows=1, usecols=range(4, colsnum)) 
    # loading targets str
    targets = jf.loadTargetsStrFromFirstCol(filename)
    # names of samples
    samples = temp[4:] 

    print 'SVD...'
    U, S, Vt = np.linalg.svd(data, full_matrices=False)
    
    index = S < 0.7 * np.mean(S)
    new_S = np.diag(S * index)
    
    # reconstruct data matrix
    data_new = np.dot(U, np.dot(new_S, Vt))

    # save to files
    file_u = open(output + '.U', 'w')
    file_s = open(output + '.S', 'w')
    file_vt = open(output + '.Vt', 'w')
    
    print 'Saving SVD files...'
    np.savetxt(file_u, U, delimiter='\t')
    np.savetxt(file_s, S, delimiter='\t')
    np.savetxt(file_vt, Vt, delimiter='\t')
    file_u.close()
    file_s.close()
    file_vt.close()

    print 'Saving matrix..'
    jf.saveRPKMMatrix(output + '.SVD', samples, targets, np.transpose(data_new))

def discover(args) :
    paramsfile = args.params
    sample_req = args.sample
    hetsnp = args.hetsnp
    tagsnp = args.tagsnp
    vcf_file = args.vcf

    if hetsnp == 'True' or hetsnp == 'TRUE':
        hetsnp = True
    else:
        hetsnp = False
    
    if tagsnp == 'True' or tagsnp == 'TRUE':
        tagsnp = True
    else:
        tagsnp = False

    datafile = args.rpkm_matrix
    f_dir = os.path.dirname(datafile)
    if f_dir != '':
        f_dir = f_dir + '/'

    if args.output:
		outputfile = f_dir + str(args.output)

    tagsnp_file = args.tagsnp_file
    mode = args.mode

    sample_flag = False #used to check whether sample_req exists

    # Build a reference set 
    if mode == 'single' or mode == 'baseline' or mode == 'reference' or mode == 'ref':
        print 'Building the reference dataset...'
        dataloader = DataManager(datafile)
        samples_np = dataloader.getAllSamples()
        dataloader.closeFile()
        print 'Baseline is Done.'

    print 'Loading data file...',
    dataloader = DataManager(datafile)
    print 'Done!'
    print 'Loading paramters...',
    params = dataloader.getParams(paramsfile)
    print 'Done!'
    dataloader.skipHeadline()
    sample = dataloader.getNextSample()

    targets_list = dataloader.getTargetsList()
    output_aux = file(outputfile+'.aux', 'w')
    output_aux.write('SAMPLE_ID\tCNV_TYPE\tFULL_INTERVAL\tINDEX\tINTERVAL\tREAD_DEPTH\n')
    output = file(outputfile,'w')
    output.write('SAMPLE_ID\tCNV_TYPE\tINTERVAL\tCHROMOSOME\tSTART\tSTOP\tLENGTH\n')

    if (hetsnp or tagsnp) and vcf_file == '':
        print 'Error: please indicate a vcf file!'
        system.exit(0)

    if vcf_file != '':
        vcf_reader = VCFReader(vcf_file)
    else:
	vcf_reader = False

    if tagsnp:
        print 'Loading tagSNP information ...',
        cnp_dict = vcf_reader.loadTagSNP(tagsnp_file)
        print 'Done!'

    while sample :
        if sample_req == '' or (sample_req != '' and sample['sample_id'] == sample_req):
            sample_flag = True
            print time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) ,sample_req,'......'

            #Renjie added: To check whether the VCF contains sample_req.
            vcf_checker = vcf.Reader(open(vcf_file,'r'))
            if sample['sample_id'] in vcf_checker.samples:
                sample_in_VCF = True
            elif sample_req in vcf_checker.samples:
                sample_in_VCF = True
            else:
                print 'No sample %s in VCF file.'%sample_req
                sample_in_VCF = False

            if hetsnp and sample_in_VCF :
                print 'Parsing SNV information from VCF file for: ' + sample['sample_id']
                snp_info = vcf_reader.getSNPInfo(sample['sample_id'], targets_list)

            if tagsnp and sample_in_VCF:
                print 'Analysing tagSNP information from tagSNP database for: ' + sample['sample_id'],
                cnp_list = vcf_reader.findTagSNPForSample(sample['sample_pop'], sample['sample_id'], cnp_dict)
                tagsnp_info_list = vcf_reader.findExonWithTagSNP(cnp_list, targets_list, overlap_threshold=0.5)
                print len(tagsnp_info_list)

            #estimate NB paramters from sample['observations']  
            sample_observations = []
            remove_list = []
            sample['observations'] = [ float(x) for x in sample['observations']]
            
            #slicing: target_index is used to split observations sequence
            target_index_begin = 0
            target_index_end = 0
            temp = 1

            sample_observations_list = []
            snp_info_list = []

            for i, targets in enumerate(targets_list):
                target_index_end = target_index_begin + len(targets)
                if hetsnp and sample_in_VCF:
                    snp_info_list.append(snp_info[target_index_begin:target_index_end])
                sample_observations_list.append(sample['observations'][target_index_begin:target_index_end])

                target_index_begin = target_index_end

            # Filtering:
            if mode == 'svd' or mode == 'SVD' or mode == 'pooled' or mode == 'pooled-sample':
                for i in range(len(sample_observations_list)):
                    sample_observations_list[i] = ndarray.tolist(stats.zscore(sample_observations_list[i]))

            elif mode == 'baseline' or mode == 'reference'  or mode == 'single' or mode == 'single-sample':
                # filtering lists whose observation equals to 0

                for i in range(len(targets_list)):
                    rem_index = []
                    for j in range(len(targets_list[i])):
                        value = sample_observations_list[i][j]
                        if np.isnan(float(value)):
                            rem_index.append(j)
                    #filter target_list, snp_list and observation_list    
                    targets_list[i] = jf.filter_list_by_list(targets_list[i], rem_index)
                    sample_observations_list[i] = jf.filter_list_by_list(sample_observations_list[i], rem_index)
                    if hetsnp and sample_in_VCF:
                        snp_info_list[i] = jf.filter_list_by_list(snp_info_list[i], rem_index)
                    if tagsnp and sample_in_VCF:
                        tagsnp_info_list[i] = jf.filter_list_by_list(tagsnp_info_list[i], rem_index)

                #Parameters estimation
                observations_all_list = []
                for i in range(len(sample_observations_list)):
                    observations_all_list.extend(sample_observations_list[i])

                parameterLoader = ParameterEstimation(observations_all_list)
                parameterList = parameterLoader.fit(observations_all_list,0.01,0.99)
                print "Estimated Paramters: ",parameterList
                params.append(parameterList[0])#mu
                params.append(parameterList[1])#sd
                    
            for i, targets in enumerate(targets_list):
                print 'Running HMM for sample[' + sample['sample_id'] + ']: ',
                print 'chr' + targets[0]._chr + ' [' + str(temp) + '|' + str(len(targets_list)) + ']'
                temp += 1
		
                #Run the HMM 
                if not hetsnp and not tagsnp:
                    modelParams = ModelParams(mode, params, targets, het_nums=0, tagsnp=0)
                elif sample_in_VCF and hetsnp and not tagsnp:
                	modelParams = ModelParams(mode, params, targets, snp_info_list[i], tagsnp=0)
                elif sample_in_VCF and not hetsnp and tagsnp:
                	modelParams = ModelParams(mode, params, targets, het_nums=0, tagsnp=tagsnp_info_list[i])
                elif sample_in_VCF and hetsnp and tagsnp:
                	modelParams = ModelParams(mode, params, targets, snp_info_list[i], tagsnp_info_list[i])
                elif not sample_in_VCF and hetsnp and tagsnp:
                    modelParams = ModelParams(mode, params, targets, het_nums=0, tagsnp=0)
                else:
                    pdb.set_trace()
	
                model = Model(mode, modelParams, sample_observations_list[i])
                pathlist = list()
                
                if vcf_reader and sample_in_VCF:
                    pathlist = model.forwardBackward_Viterbi(mode, if_snp = True)
                else:
                    pathlist = model.forwardBackward_Viterbi(mode, if_snp = False)
                dataloader.outputCNVaux(output_aux, sample['sample_id'], targets, pathlist, sample_observations_list[i])
                dataloader.outputCNV(output, sample['sample_id'], targets, pathlist, sample_observations_list[i])

        sample = dataloader.getNextSample()

    output.close()
    output_aux.close()
    dataloader.closeFile()

    if not sample_flag:
        print 'Could not find the sample_id specified.'

def merge_results(args):
    datafile_svd = args.datafile_svd
    datafile_dis = args.datafile_dis
    output = args.output
    conflict_output = output + '.conflict_results'
    print conflict_output
    try:
        svd_file = open(datafile_svd, 'r')
        dis_file = open(datafile_dis, 'r')
        output_file = open(output, 'w')
        conflict_file = open(conflict_output, 'w')
    except IOError:
        sys.exit(0)

    conflict_file.write('-'*80 + '\n')
    # skip first row
    output_file.write(svd_file.readline().strip('\n') + '\tmode' + '\n')
    dis_file.readline()

    svd_line = svd_file.readline()
    dis_line = dis_file.readline()

    while svd_line != '' or dis_line != '':
        if svd_line == '' and dis_line != '':
            output_file.write(dis_line.strip('\n') + '\tdistribution\n')
            dis_line = dis_file.readline()
            continue
        if svd_line != '' and dis_line == '':
            output_file.write(svd_line.strip('\n') + '\tsvd\n')
            svd_line = svd_file.readline()
            continue

        svd_temp = svd_line.strip('\n').split('\t')
        dis_temp = dis_line.strip('\n').split('\t')
        
        if svd_temp[3] == 'X':
            svd_chr = 23
        elif svd_temp[3] == 'Y':
            svd_chr = 24
        else:
            svd_chr = int(svd_temp[3])
        svd_start = int(svd_temp[4])
        svd_stop = int(svd_temp[5])

        if dis_temp[3] == 'X':
            dis_chr = 23
        elif dis_temp[3] == 'Y':
            dis_chr = 24
        else:
            dis_chr = int(dis_temp[3])
        dis_start = int(dis_temp[4])
        dis_stop = int(dis_temp[5])

        if svd_chr < dis_chr:
            output_file.write(svd_line.strip('\n') + '\tsvd\n')
            svd_line = svd_file.readline()
            continue
        elif svd_chr > dis_chr:
            output_file.write(dis_line.strip('\n') + '\tdistribution\n')
            dis_line = dis_file.readline()
            continue
        else:
            if svd_stop < dis_start:
                output_file.write(svd_line.strip('\n') + '\tsvd\n')
                svd_line = svd_file.readline()
                continue 
            elif svd_start > dis_stop:
                output_file.write(dis_line.strip('\n') + '\tdistribution\n')
                dis_line = dis_file.readline()
                continue	

            else:
                temp = [''] * 8
                temp[0] = svd_temp[0]
                if svd_temp[1] != dis_temp[1]:
                    print "Warning: conflict exists between two methods"
                    conflict_file.write(svd_line.strip('\n') + '\tsvd\n')
                    conflict_file.write(dis_line.strip('\n') + '\tdistribution\n')
                    conflict_file.write('-'*80 + '\n')
                    svd_line = svd_file.readline()
                    dis_line = dis_file.readline()
                    continue
 
                else:
                    temp[1] = svd_temp[1]
                start = svd_start if svd_start <= dis_start else dis_start
                stop = svd_stop if svd_stop >= dis_stop else dis_stop
                if svd_chr == 23:
                    svd_chr = 'X'
                elif svd_chr == 24:
                    svd_chr = 'Y'
                temp[2] = str(svd_chr) + ':' + str(start) + ':' + str(stop)
                temp[3] = str(svd_chr)
                temp[4] = str(start)
                temp[5] = str(stop)
                temp[6] = str(stop - start + 1)
                temp[7] = 'svd[' + svd_temp[2] + '] | distribution[' + dis_temp[2] + ']'
                output_file.write('\t'.join(temp) + '\n')
                svd_line = svd_file.readline()
                dis_line = dis_file.readline()

    output_file.close()
    conflict_file.close()


parser = argparse.ArgumentParser(prog='jxcnv', description='Designed by jx.')
subparsers = parser.add_subparsers()

#BAM List -> RPKM
svd_parser = subparsers.add_parser('rpkm', help="Create RPKM matrix from a BAM list")
svd_parser.add_argument('--target', required=True, help='Target definition file')
svd_parser.add_argument('--maq', required=False, type=int, default=20, help='MAQ threshold')
svd_parser.add_argument('--input', required=True, help='BAM file list, each line for each sample')
svd_parser.add_argument('--output', required=True, help='Directory for RPKM files')
svd_parser.set_defaults(func=bamlist2RPKM)

#normalize read depth signal
svd_parser = subparsers.add_parser('normalize_gc', help="GC correction for a given sample")
svd_parser.add_argument('--windows', required=True, help='Please input the target information including GC content')
svd_parser.add_argument('--input', required=True, help='Please input a read depth file for a given sample')
svd_parser.add_argument('--output', required=False, help='GC corrected RD signal')
svd_parser.set_defaults(func=normalize_gc)

#RPKM files -> Matrix
svd_parser = subparsers.add_parser('merge_rpkm', help="Merge RPKM files to a matrix")
svd_parser.add_argument('--rpkm_dir', required=True, help='RPKM files')
svd_parser.add_argument('--target', required=True, help='Target definition file')
svd_parser.add_argument('--output', required=False, help='Matrix file')
svd_parser.set_defaults(func=RPKM2Matrix)

# Filter matrix by GC content, mapping ability and exon length
svd_parser = subparsers.add_parser('filter', help="Filter matrix by GC content, mapping ability and exon length")
svd_parser.add_argument('--rpkm_matrix', required=True, help='Matrix of RPKM values')
svd_parser.add_argument('--ref_file', required=False, help='Reference file for the calculation of GC percentage')
svd_parser.add_argument('--map_file', required=False, help='Mapping ability file.')
svd_parser.add_argument('--filter_params', required=True, help='Parameters of filtering')
svd_parser.add_argument('--output', required=False, help='Filtered matrix')
svd_parser.set_defaults(func=filter_rpkm)

#SVD
svd_parser = subparsers.add_parser('svd', help="SVD")
svd_parser.add_argument('--rpkm_matrix', required=True, help='')
svd_parser.add_argument('--output', required=False, help='')
# svd_parser.add_argument('--svd', type=int, required=True, help='Number of components to remove')
svd_parser.set_defaults(func=svd)

#CNV discover
cnv_parser = subparsers.add_parser('discover', help="Run HMM to discover CNVs")
cnv_parser.add_argument('--params', required=True, help='Parameters used by HMM')
cnv_parser.add_argument('--rpkm_matrix', required=True, help='RPKM matrix.')
cnv_parser.add_argument('--mode',required=True, default='SVD', help='Data normalization by SVD or baseline mode.')
cnv_parser.add_argument('--output', required=True, help='Output file.')
cnv_parser.add_argument('--sample', required=False, default='', help='Optionally, users can choose one sample to run.')
cnv_parser.add_argument('--vcf', required=False, default='', help='Optionally, users can input snp information by specifing a vcf file')
cnv_parser.add_argument('--hetsnp', required=False, default=False)
#cnv_parser.add_argument('--no-hetsnp', dest='hetsnp', action='store_true')
#cnv_parser.set_defaults(hetsnp=True)
cnv_parser.add_argument('--tagsnp', required=False, default=False)
#cnv_parser.add_argument('--no-tagsnp', dest='tagsnp', action='store_true')
#cnv_parser.set_defaults(tagsnp=True)
cnv_parser.add_argument('--tagsnp_file',required = False, help='TagSNP file location.')
cnv_parser.set_defaults(func=discover)

#Merge results
cnv_parser = subparsers.add_parser('merge', help="Merge results from different methods")
cnv_parser.add_argument('--datafile_svd', required=True)
cnv_parser.add_argument('--datafile_dis', required=True)
cnv_parser.add_argument('--output', required=True)
cnv_parser.set_defaults(func=merge_results)

args = parser.parse_args()
args.func(args)
