#!/usr/bin/env python
# coding: utf-8

# # CNN Prediction

from __future__ import print_function
import os
import re
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

import function_dl as func_dl


img_width, img_height = 224, 224


## For pcgc experimental cnvs
cnv_info_file = '/home/rt2776/cnv_espresso/predict_pcgc/images_new/pcgc_NimbleGenV2_data_withImagePath.csv'
cnv_info_df   = pd.read_csv(cnv_info_file)

entire_cnv_images_path_list  = cnv_info_df['entire_cnv_path']

# ### Loading images from list to numpy array
img_np = func_dl.loadImgs(entire_cnv_images_path_list, img_width, img_height)
pdb.set_trace()
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
model_path = '/home/rt2776/cnv_espresso/images_rare_3classes/data_backup/model_h5/rare_entire_cnv_MobileNet_v1_fine-tuning_3classes.h5'
print("Loading %s ... from %s"%(model_name, model_path))
MobileNet_model = keras.models.load_model(model_path, custom_objects=custom_objects)

img_pred = MobileNet_model.predict(img_np)

pred_output_df = copy.deepcopy(cnv_info_df)
pred_output_df.shape

pred_output_df.insert(pred_output_df.shape[1], 'Prob_True', "")
pred_output_df.insert(pred_output_df.shape[1], 'Prob_False', "")
pred_output_df.insert(pred_output_df.shape[1], 'Prediction', "")
pred_output_df.insert(pred_output_df.shape[1], 'Pred_status', "")

## output to file
output_path = '/home/rt2776/cnv_espresso/predict_pcgc/'
pred_output_df.to_csv(output_path+'pcgc_NimbleGenV2_data_prediction_0223.csv',index=False)

