from __future__ import print_function
import os
import re
import copy
import random
import datetime
import timeit
import pdb
import PIL
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import classification_report
from sklearn import metrics
from sklearn.model_selection import KFold, StratifiedKFold
import sklearn

import tensorflow as tf
from tensorflow import keras
import keras.preprocessing
from keras.models import Sequential, Model
#from keras.utils import to_categorical #tf2.5.0 does not support
from keras.layers import Dense, Conv2D, MaxPooling2D, Dropout, Flatten
from keras.models import load_model
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint
from keras import backend

import function_dl as func_dl
import function as func


def cnn_train(true_del_file, true_dup_file, false_del_file, false_dup_file, use_gpu, batch_size, epochs, output_model_dir):
    # GPU or CPU selection
    if use_gpu == False or use_gpu == 'False':
        print("Using CPU ...")
    else:
        os.environ['CUDA_VISIBLE_DEVICES'] = "1"
        physical_devices = tf.config.experimental.list_physical_devices('GPU') 
        physical_devices

        config = tf.compat.v1.ConfigProto()
        config.gpu_options.allow_growth = True

    img_width, img_height = 224, 224
    seed = 2021

    # check output_model_dir
    os.makedirs(output_model_dir, exist_ok=True)

    # ## Importing data from scratch
    ## For rare CNVs

    true_del_df  = pd.read_csv(true_del_file,  header=0,sep='\t')
    false_del_df = pd.read_csv(false_del_file, header=0,sep='\t')

    true_dup_df  = pd.read_csv(true_dup_file,  header=0,sep='\t')
    false_dup_df = pd.read_csv(false_dup_file, header=0,sep='\t')

    true_del_images_path_list  = true_del_df['img_path']
    false_del_images_path_list = false_del_df['img_path']

    true_dup_images_path_list  = true_dup_df['img_path']
    false_dup_images_path_list = false_dup_df['img_path']

    print("The shape of each type:")
    print("True  DEL:", true_del_images_path_list.shape)
    print("True  DUP:", true_dup_images_path_list.shape)
    print("False DEL:", false_del_images_path_list.shape)
    print("False DUP:", false_dup_images_path_list.shape)


    # ### Loading images from list to numpy array
    func.showDateTime(end='\t')
    print("Loading images from list to numpy array ...")
    true_del_img_np = func_dl.loadImgs(true_del_images_path_list, img_width, img_height)
    true_del_img_np.shape

    false_del_img_np = func_dl.loadImgs(false_del_images_path_list, img_width, img_height)
    false_del_img_np.shape

    true_dup_img_np = func_dl.loadImgs(true_dup_images_path_list, img_width, img_height)
    true_dup_img_np.shape

    false_dup_img_np = func_dl.loadImgs(false_dup_images_path_list, img_width, img_height)
    false_dup_img_np.shape


    # ### Generate labels for entire CNVs
    # Three classes
    true_del_label  = [0 for i in range(0,len(true_del_img_np))]
    false_del_label = [1 for i in range(0,len(false_del_img_np))]
    true_dup_label  = [2 for i in range(0,len(true_dup_img_np))]
    false_dup_label = [1 for i in range(0,len(false_dup_img_np))]

    print(len(true_del_label), len(false_del_label), len(true_dup_label), len(false_dup_label))

    # ### Combine true & false data for entire CNVs
    combined_cnv_info_df = true_del_df.append(false_del_df, ignore_index=True)
    combined_cnv_info_df = combined_cnv_info_df.append(true_dup_df, ignore_index=True)
    combined_cnv_info_df = combined_cnv_info_df.append(false_dup_df, ignore_index=True)

    combined_img = np.vstack((true_del_img_np, false_del_img_np, true_dup_img_np, false_dup_img_np))

    combined_label = true_del_label + false_del_label + true_dup_label + false_dup_label
    len(combined_label)

    # ## Normalization
    # Find the shape of input images and create the variable input_shape
    nRows,nCols,nDims = combined_img.shape[1:]
    input_shape = (nRows, nCols, nDims)
    print("The shape of input tensor:",input_shape)

    # Change to float datatype
    combined_img = combined_img.astype('float32')

    # Scale the data to lie between 0 to 1
    combined_img /= 255

    # Change the labels from integer to categorical data
    combined_label_one_hot = tf.keras.utils.to_categorical(combined_label)


    # ## Find the unique numbers from the train labels
    classes = np.unique(combined_label)
    nClasses = len(classes)
    print('Total number of outputs : ', nClasses)
    print('Output classes : ', classes)
    print("3 classes label: 0-True deletion; 1-Diploid (False del & False dup); 2-True duplication")


    # ## Train the deep nerual model by hold-out validation
    func.showDateTime(end='\t')
    print("Split dataset into training (60%), validation (20%) and testing (20%) dataset ...")

    ## split image arrays
    train_img, test_img, train_label, test_label, train_cnv_info_df, test_cnv_info_df = train_test_split(combined_img,
                                                                                                        combined_label_one_hot,
                                                                                                        combined_cnv_info_df,
                                                                                                        test_size=0.2,
                                                                                                        shuffle=True,
                                                                                                        random_state=seed)

    train_img, val_img, train_label, val_label, train_cnv_info_df, val_cnv_info_df = train_test_split(train_img,
                                                                                                      train_label,
                                                                                                      train_cnv_info_df,
                                                                                                      test_size=0.25,
                                                                                                      shuffle=True,
                                                                                                      random_state=seed) # 0.25*0.8=0.2

    combined_img.shape, train_img.shape, val_img.shape, test_img.shape
    combined_label_one_hot.shape, train_label.shape, val_label.shape, test_label.shape


    # ## CNN (Transfer learning and fine-tuning)
    # ### Using the pretrained MobileNet v1 architecture
    # - Firstly, we keep all the weights of base model frozen to train the FC layers.
    func.showDateTime(end='\t')
    print("CNN transfer learning and fine-tuning ...")
    print("\t 1. frozen the pretrained MobileNet v1 and train the fully connected layers ...")

    model_name='MobileNet_v1_fine_tuning'
    base_model = tf.keras.applications.MobileNet(
                                    weights='imagenet', # Load weights pre-trained model.
                                    input_shape=(224, 224, 3),  
                                    include_top=False)  # Do not include the ImageNet classifier at the top.

    base_model.trainable = False
    inputs = keras.Input(shape=(224, 224, 3)) 
    x = base_model(inputs, training=False)

    # Convert features of shape `base_model.output_shape[1:]` to vectors
    x = keras.layers.GlobalAveragePooling2D()(x)
    # A Dense classifier with a single unit (binary classification)
    outputs = keras.layers.Dense(nClasses,activation='softmax')(x)
    model   = keras.Model(inputs, outputs)
    model.summary()

    model.compile(optimizer=keras.optimizers.Adam(),
                  loss='categorical_crossentropy',
                  metrics=['accuracy', func_dl.f1_m, func_dl.precision_m, func_dl.recall_m])

    print("Training by MobileNet_v1 model ...")

    model_file = output_model_dir + "/" +  model_name + "_" + str(nClasses) + "classes.h5"

    es = EarlyStopping(monitor  ='val_loss', mode='min', verbose=1, patience=3)
    mc = ModelCheckpoint(model_file,
                         monitor='val_accuracy',
                         mode   ='max', 
                         verbose=1, 
                         save_best_only=True)

    history = model.fit(train_img, train_label,
                        batch_size = batch_size, 
                        epochs =epochs,
                        verbose=1, 
                        validation_data=(val_img, val_label), 
                        callbacks=[es, mc])

    print("\n")
    loss, accuracy, f1_score, precision, recall = model.evaluate(test_img, test_label)

    func_dl.draw_loss_accuracy_curves(history, "CNN transfer learning", 
                                                output_img_file=output_model_dir + "/"+ "loss_accuracy_curves_1.pdf")
    func_dl.confusion_matrix(model, test_img, test_label, nClasses, 
                                                output_img_file=output_model_dir + "/"+ "confusion_matrix_1.pdf")
    fpr, tpr, thresholds, auc = func_dl.pred_roc_data(model, test_img, test_label) 
    func_dl.draw_single_roc_curve(tpr, fpr, auc, output_img_file=output_model_dir + "/"+ "ROC_curve_1.pdf")


    # ### Fine-tuning
    # - Secondly, Once your model has converged on our train data, we unfreeze all or part of the base model and retrain the whole model end-to-end with a very low learning rate.

    func.showDateTime(end='\t')
    print("\t 2. Fine tuning the MobileNet_v1 model ...")
    model_file = output_model_dir + "/" + model_name + "_" + str(nClasses) + "classes.h5"

    base_model.trainable=True
    model.summary()

    model.compile(optimizer=keras.optimizers.Adam(1e-5),
        loss='categorical_crossentropy', metrics=['accuracy', func_dl.f1_m, func_dl.precision_m, func_dl.recall_m])

    es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=3)
    mc = ModelCheckpoint(model_file,
                         monitor='val_accuracy',
                         mode   ='max', 
                         verbose=1, 
                         save_best_only=True)

    history = model.fit(train_img, train_label,
                        batch_size = batch_size, 
                        epochs  = epochs,
                        verbose = 1, 
                        validation_data = (val_img, val_label), 
                        callbacks = [es, mc])
    print("\n")
    loss, accuracy, f1_score, precision, recall = model.evaluate(test_img, test_label)

    func_dl.draw_loss_accuracy_curves(history, "CNN fine-tuning", 
                                                output_img_file=output_model_dir + "/"+ "loss_accuracy_curves_2.pdf")
    func_dl.confusion_matrix(model, test_img, test_label, nClasses, 
                                                output_img_file=output_model_dir + "/"+ "confusion_matrix_2.pdf")
    fpr, tpr, thresholds, auc = func_dl.pred_roc_data(model, test_img, test_label) 
    func_dl.draw_single_roc_curve(tpr, fpr, auc, output_img_file=output_model_dir + "/"+ "ROC_curve_2.pdf")

    func.showDateTime()
    print("[Done]. Please check the trained model at",model_file)

