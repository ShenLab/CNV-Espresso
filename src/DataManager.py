from hmm.Interval import *
import os
import numpy as np
import vcf
import pdb

class DataManager(object) :
    def __init__(self, filename) :
        self._filename = None 
        if os.path.exists(filename) :
            self._filename = filename
            self._datafile = open(self._filename)
            self._hasSkipHeadline = False
        else:
            raise Exception('Data file does not exist')
        
    def getTargetsList(self) :
        if self._filename :
            datafile = open(self._filename)
            headline = datafile.readline().strip()
            targets = []
            #since there could be more than one chrs in a target line
            #the targets of different chrs should be put in different list
            targets_list = []
            headline_temp = headline.split('\t')[1:]
            pre_str = headline_temp[0]
            for targetstr in headline_temp:
                if targetstr.split(':')[0] == pre_str.split(':')[0]:
                    targets.append(self.buildIntervalFromStr(targetstr))
                else:
                    targets_list.append(targets)
                    targets = [self.buildIntervalFromStr(targetstr)]
                pre_str = targetstr
            targets_list.append(targets)

            return targets_list

    def buildIntervalFromStr(self, s):
        temp1 = s.split(':')
        _chr = temp1[0]
        temp2 = temp1[1].split('-')
        bp1 = temp2[0]
        bp2 = temp2[1]
        interval = Interval(_chr, bp1, bp2)
        
        return interval
        
    def skipHeadline(self):
        if self._datafile and not self._hasSkipHeadline :
            self._datafile.readline()
            self._hasSkipHeadline = True

    def getNextSample(self):
        samplestr = self._datafile.readline()
        if samplestr :
            samplestr.strip('\n')
            samlist = samplestr.split('\t')
            sample = {'sample_id':samlist[0].split('.')[0], 'sample_pop':samlist[0].split('.')[1], 'observations':samlist[1:]}
            
            return sample
        return None

    def getSampleNum(self):
        sample_num = 0
        sample = self._datafile.readline()
        while sample:
            sample_num += 1
            sample = self._datafile.readline()

        return sample_num

    def getAllSamples(self):
        # skip 1st row and 1st column
        line = self._datafile.readline()
        line = line.strip()
        colsnum = len(line.split('\t'))

        data = np.loadtxt(self._datafile, dtype=np.float, delimiter='\t', skiprows=0, usecols=range(1, colsnum)) 
        return data

    def getParams(self, paramsfile):
        try:
            pfile = open(paramsfile)
            pline = pfile.readline().strip('\n')
            params = pline.split('\t')
            pfile.close()
            
            return params
        except IOError, msg:
            print msg
            return None
        
    def closeFile(self):
        if self._datafile:
            self._datafile.close()

    def outputCNVaux(self, output, sample_id, targets, pathlist, observations) :
        dup_begin_index = 0
        dup_begin_flag = False

        del_begin_index = 0
        del_begin_flag = False
        
        for i in range(len(targets)) :
            if dup_begin_flag and (pathlist[i] != 'DUP' or i+1 == len(targets)):
                full_interval = targets[dup_begin_index]._chr + ':' + str(targets[dup_begin_index]._bp1) + '-' + str(targets[i-1]._bp2)
                for j in range(dup_begin_index-2 if dup_begin_index-2>0 else 0, dup_begin_index):
                    output.write(sample_id + '\t' + pathlist[dup_begin_index] + '\t' + full_interval + '\tU-' + str(dup_begin_index-j) + '\t' + targets[j].getInfo() + '\t' + str(observations[j]) +'\n')
                for j in range(dup_begin_index, i):
                    output.write(sample_id + '\t' + pathlist[j] + '\t' + full_interval + '\t' + str(j-dup_begin_index+1) + '\t' + targets[j].getInfo() + '\t' + str(observations[j]) + '\n')
                for j in range(i, i+2 if i+2<len(targets) else len(targets)):
                    output.write(sample_id + '\t' + pathlist[dup_begin_index] + '\t' + full_interval + '\tD+' + str(j-i+1) + '\t' + targets[j].getInfo() + '\t' + str(observations[j]) + '\n')
                dup_begin_flag = False

            if del_begin_flag and (pathlist[i] != 'DEL' or i+1 == len(targets)):
                full_interval = targets[del_begin_index]._chr + ':' + str(targets[del_begin_index]._bp1) + '-' + str(targets[i-1]._bp2)
                for j in range(del_begin_index-2 if del_begin_index-2>0 else 0, del_begin_index):
                    output.write(sample_id + '\t' + pathlist[del_begin_index] + '\t' + full_interval + '\tU-' + str(del_begin_index-j) + '\t' + targets[j].getInfo() + '\t' + str(observations[j]) +'\n')
                for j in range(del_begin_index, i):
                    output.write(sample_id + '\t' + pathlist[j] + '\t' + full_interval + '\t' + str(j-del_begin_index+1) + '\t' + targets[j].getInfo() + '\t' + str(observations[j])  + '\n')
                for j in range(i, i+2 if i+2<len(targets) else len(targets)):
                    output.write(sample_id + '\t' + pathlist[del_begin_index] + '\t' + full_interval + '\tD+' + str(j-i+1) + '\t' + targets[j].getInfo() + '\t' + str(observations[j]) + '\n')
                del_begin_flag = False

            if pathlist[i] == 'DUP' and not dup_begin_flag:
                dup_begin_index = i
                dup_begin_flag = True

            if pathlist[i] == 'DEL' and not del_begin_flag:
                del_begin_index = i
                del_begin_flag = True 

    def outputCNV(self, output, sample_id, targets, pathlist, observations) :
        dup_begin_index = 0
        dup_begin_flag = False

        del_begin_index = 0
        del_begin_flag = False
        
        for i in range(len(targets)) :
            if dup_begin_flag and pathlist[i] != 'DUP':
                full_interval = targets[dup_begin_index]._chr + ':' + str(targets[dup_begin_index]._bp1) + '-' + str(targets[i-1]._bp2)

                dup_begin_flag = False
                try:
                    cnv_lenth = targets[i-1]._bp2 - targets[dup_begin_index]._bp1 + 1
                    output.write(sample_id + '\tDUP\t' + full_interval + '\t' + targets[dup_begin_index]._chr + '\t' + str(targets[dup_begin_index]._bp1) + '\t' + str(targets[i-1]._bp2) +'\t'+ str(cnv_lenth) + '\n')
                except Exception, e:
                    pass

            if del_begin_flag and pathlist[i] != 'DEL':
                full_interval = targets[del_begin_index]._chr + ':' + str(targets[del_begin_index]._bp1) + '-' + str(targets[i-1]._bp2)

                del_begin_flag = False
                try:
                    cnv_lenth = targets[i-1]._bp2 - targets[del_begin_index]._bp1 + 1
                    output.write(sample_id + '\tDEL\t' + full_interval + '\t' + targets[del_begin_index]._chr + '\t' + str(targets[del_begin_index]._bp1) + '\t' + str(targets[i-1]._bp2) +'\t'+ str(cnv_lenth) + '\n')
                except Exception, e:
                    pass

            if pathlist[i] == 'DUP' and not dup_begin_flag:
                dup_begin_index = i
                dup_begin_flag = True

            if pathlist[i] == 'DEL' and not del_begin_flag:
                del_begin_index = i
                del_begin_flag = True 


                    







