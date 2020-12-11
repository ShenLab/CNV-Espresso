import numpy as np
import pysam
import pdb
import bx.bbi.bigwig_file
import sys
import os

def loadMatrixFromFile(filename, skiprows, skipcols, type='float', delimiter='\t'):
    return

def loadTargets(target_filename):
	#Renjie modified

    targetfile = open(target_filename)
    targets = []
    target_id = 1

    for line in targetfile.readlines():
        line = line.strip('\n')
        temp = line.split('\t')
		
        targets.append({'chr': temp[0], 'start': int(temp[1]), 'stop': int(temp[2])})

    return targets

def loadTargetsStr(target_filename):
    targetfile = open(target_filename)
    targets = []

    for line in targetfile.readlines():
        line = line.strip()
        temp = line.split('\t')
        targets.append(temp[0]+':'+temp[1]+'-'+temp[2])

    return targets

def loadTargetsFromFirstCol(filename):
    f = open(filename)
    targets_str = np.loadtxt(f, dtype=np.str, delimiter='\t', skiprows=1, usecols=(0,)) 
    targets = []
    for line in targets_str:
        chr = line.split('\t')[0].split(':')[0]
        start = int(line.split('\t')[0].split(':')[1].split('-')[0])
        stop = int(line.split('\t')[0].split(':')[1].split('-')[1])
        targets.append({'chr': chr, 'start': start, 'stop': stop})

    return {'targets': targets, 'targets_str': targets_str}

def loadTargetsStrFromFirstCol(filename):
    f = open(filename)
    return np.loadtxt(f, dtype=np.str, delimiter='\t', skiprows=1, usecols=(0,)) 

def loadValuesByCol(filename, col):
    if not os.path.exists(filename):
        print filename, 'does not exists!'
        sys.exit(0)
    return np.loadtxt(open(filename), dtype=np.int, delimiter='\t', skiprows=1, usecols=(col,))

def loadRPKMMatrix(filename):
    f = open(filename)
    line = f.readline()
    line = line.strip()
    colsnum = len(line.split('\t'))
    samples = line.split('\t')[4:]
    
    f.seek(0)
    rpkm = np.loadtxt(f, dtype=np.float, delimiter='\t', skiprows=1, usecols=range(4, colsnum)) 
    f.seek(0)
    annotation = np.loadtxt(f, dtype=np.int, delimiter='\t', skiprows=1, usecols=range(1, 4)) 
    f.seek(0)
    targets_str = np.loadtxt(f, dtype=np.str, delimiter='\t', skiprows=1, usecols=(0,)) 
    targets = []
    for line in targets_str:
        chr = line.split('\t')[0].split(':')[0]
        start = int(line.split('\t')[0].split(':')[1].split('-')[0])
        stop = int(line.split('\t')[0].split(':')[1].split('-')[1])
        targets.append({'chr': chr, 'start': start, 'stop': stop})
    
    return {'samples': samples, 'rpkm': rpkm, 'annotation': annotation, 'targets': targets, 'targets_str': targets_str}

def saveRPKMMatrix(filename, samples, targets, rpkm):
    rpkm_f = open(filename, 'w')
    rpkm_f.write('Matrix\t' + '\t'.join(targets) + '\n')
    if len(samples) != len(rpkm):
        print 'Error when saving rpkm. The row numbers do not match'
        sys.exit(0)

    for i in range(len(samples)):
        rpkm_f.write(samples[i] + '\t' + '\t'.join(str(r) for r in rpkm[i]) + '\n')
    rpkm_f.close()

def loadNormValues(filename):
    #Renjie revised
    return np.loadtxt(open(filename), dtype=np.int, delimiter='\t', skiprows=1, usecols=(1,)) 

def saveNormValues(filename, targets, values, header):
    f = open(filename, 'w')
    f.write('Targets\t' + header + '\n')
    if len(targets) != len(values):
        print 'Error when saving normalization values. The row numbers do not match'
        sys.exit(0)

    for i in range(len(targets)):
        f.write(targets[i] + '\t' + str(values[i]) + '\n')
    f.close()
    

def chrInt2Str(chromosome_int):
    return 'chr' + str(chromosome_int)
    if int(chromosome_int) == 23:
        return 'chrX'
    elif int(chromosome_int) == 24:
        return 'chrY'
    else:
        return 'chr' + str(chromosome_int)

def calGCPercentage(targets, ref_file):
    fasta_f = pysam.FastaFile(ref_file)
    GC_percentage = []
    chr = targets[0]['chr']
    print 'Calculating GC content for targets in chr ' + chr
    
    for i in range(len(targets)):
        if(targets[i]['chr'] != chr):
            chr = targets[i]['chr']
            print 'Calculating GC content for targets in chr ' + chr
            
        r_region = fasta_f.fetch(targets[i]['chr'], targets[i]['start'], targets[i]['stop'])
        reg_len = targets[i]['stop'] - targets[i]['start'] + 1
        GC = 0
        n_num = 0

        for b in range(len(r_region)):
            if str(r_region[b]).upper() == 'G' or str(r_region[b]).upper() == 'C':
                GC += 1
            elif str(r_region[b]).upper() == 'N':
                n_num += 1

        if n_num / reg_len > 0.2:
            GC_p = -1
        else:
            GC_p = GC * 100 / reg_len

        GC_percentage.append(GC_p)
    
    return GC_percentage

def calMapAbility(targets, map_file):
    bw = bx.bbi.bigwig_file.BigWigFile(open(map_file, "rb"))
    map_ability = []
    chr = targets[0]['chr']
    print 'Calculating mapping ability for targets in chr ' + chr
    
    for i in range(len(targets)):
        t = targets[i] 

        if(t['chr'] != chr):
            chr = t['chr']
            print 'Calculate mapping ability for targets in chr ' + chr

        t = targets[i] 
        map_summary = bw.query(chrInt2Str(t['chr']), t['start'], t['stop'], 1)
        try:
            _map = int(map_summary[0]['mean']*100+0.5)
        except:
            _map = -1

        map_ability.append(_map)
    return map_ability

def calExonLength(targets):
    exon_length = []
    chr = targets[0]['chr']
    print 'Calculating exon length for targets in chr ' + chr

    for i in range(len(targets)):
        t = targets[i]

        if(t['chr'] != chr):
            chr = t['chr']
            print 'Calculating exon length for targets in chr ' + chr

        l = np.round((t['stop']-t['start']+1)/50)*50
        exon_length.append(l)
    return exon_length


def groupBy(data):
    result = {}
    exclude = []
    for ind, val in data:
        if val == -1:
            exlude.append(ind)
        elif result.has_key(val):
            result[val].apend(ind)
        else:
            result[val] = [ind]
    return {'exlude': exclude, 'dict': result}

def filter_list_by_list(list, filter):
    return [list[x] for x in range(len(list)) if x not in filter]

def output_to_file(results,file_name):
    path = os.getcwd() 

    if not os.path.exists(path):
        os.mkdir(path)

    file_name = path + '/' + file_name
    fp = open(file_name,'w')
    if type(results) == list:
        for result_line in results:
            if type(result_line) != list:
                fp.write(str(result_line))
            else:
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
    print 'The variable has already output to %s\n\n'%file_name
    fp.close()
