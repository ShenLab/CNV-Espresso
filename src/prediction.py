#!/usr/bin/env python
# coding: utf-8

## CNN Prediction
from __future__ import print_function
import os
import re
import sys
import pdb
import copy
import random
import datetime
import tensorflow as tf
print("Tensorflow version " + tf.__version__)
import pandas as pd
import PIL
import numpy as np
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
import sklearn
import keras
import keras.preprocessing
from keras.models import Sequential
from keras.utils import to_categorical
from keras.layers import Dense, Conv2D, MaxPooling2D, Dropout, Flatten
from keras.models import load_model
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint

import function as func
import function_dl as func_dl

# GPU selection
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
physical_devices = tf.config.experimental.list_physical_devices('GPU') 

# Variables
cnv_info_file = sys.argv[1]
model_file    = sys.argv[2] 
output_file   = sys.argv[3]

img_width, img_height = 224, 224

# Loading CNV_info and images. 
cnv_info_df   = pd.read_csv(cnv_info_file)
entire_cnv_images_path_list  = cnv_info_df['entire_cnv_path']

img_np = func_dl.loadImgs(entire_cnv_images_path_list, img_width, img_height)

# ## Normalization
# Find the shape of input images and create the variable input_shape
nRows,nCols,nDims = img_np.shape[1:]
input_shape = (nRows, nCols, nDims)
print("The shape of input tensor:",input_shape)


# Change to float datatype
img_np = img_np.astype('float32')

# Scale the data to lie between 0 to 1
img_np /= 255
# Find the unique numbers from the train labels
nClasses = 3
print('Total number classes: ', nClasses)

# ## Precision

## Load pre-calculated model
custom_objects = {"f1_m":func_dl.f1_m, "precision_m":func_dl.precision_m, "recall_m":func_dl.recall_m}

model_name = 'MobileNet_v1'
#model_path = '/home/rt2776/cnv_espresso/images_rare_3classes/data_backup/model_h5/rare_entire_cnv_MobileNet_v1_3classes.h5'
#model_path = '/home/rt2776/cnv_espresso/images_rare_3classes/data_backup/model_h5/rare_entire_cnv_MobileNet_v1_fine-tuning_3classes.h5'
print("Loading %s ... from %s"%(model_name, model_file))
MobileNet_model = keras.models.load_model(model_file, custom_objects=custom_objects)

img_pred = MobileNet_model.predict(img_np)
pred_output_df = copy.deepcopy(cnv_info_df)
pred_output_df.insert(pred_output_df.shape[1], 'Prob_DEL', "-")
pred_output_df.insert(pred_output_df.shape[1], 'Prob_DIP', "-")
pred_output_df.insert(pred_output_df.shape[1], 'Prob_DUP', "-")
pred_output_df.insert(pred_output_df.shape[1], 'Prob', "-")
pred_output_df.insert(pred_output_df.shape[1], 'Prediction', "-")
pred_output_df.insert(pred_output_df.shape[1], 'Status', "-")

num, correct_count = 0, 0
for i in range(len(img_pred)):
    num += 1
    pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prob_DEL')] = img_pred[i][0]
    pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prob_DIP')] = img_pred[i][1]
    pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prob_DUP')] = img_pred[i][2]
    pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prob')] = np.max(img_pred[i])
    
    if(np.argmax(img_pred[i]) == 0):
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] = "DEL"
    elif(np.argmax(img_pred[i]) == 1):
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] = "DIP"
    elif(np.argmax(img_pred[i]) == 2):
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] = "DUP"
    else:
        pdb.set_trace()
        
    if pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] == pred_output_df.iloc[i,pred_output_df.columns.get_loc('TYPE')]:
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Status')] = "Positive"
    else:
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Status')] = "Negative"

## output to file
pred_output_df.to_csv(output_file,index=False)
