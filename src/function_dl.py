from __future__ import print_function

import os
import re
import copy
import random
import numbers
import six
import datetime
import tensorflow as tf
print("Tensorflow version " + tf.__version__)
import pandas as pd
import pdb
import PIL
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.collections
import seaborn as sns

import sklearn
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import classification_report
from sklearn import metrics
from sklearn.model_selection import KFold, StratifiedKFold

#import keras
from tensorflow import keras
import keras.preprocessing
from keras.models import Sequential
from keras.models import Model
#from keras.utils import to_categorical
from keras.layers import Dense, Conv2D, MaxPooling2D, Dropout, Flatten
from keras.models import load_model
from keras.callbacks import EarlyStopping
from keras.callbacks import ModelCheckpoint
from keras import backend

import function as func

'''
Processing images
'''

def test(x):
    print("Testing Your input is %s"%str(x))
    
def showImg(img_data, output_image_file=None):
    if type(img_data) is list:
        image = tf.keras.preprocessing.image.load_img(img_data)
        plt.imshow(image)
    if type(img_data) is np.ndarray:
        image = tf.keras.preprocessing.image.array_to_img(img_data)
        plt.imshow(image)
        
    if output_image_file != None:
        plt.savefig(output_image_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("Output image to:",output_image_file)
        plt.close()
    
def resizeCropImg(img_file, target_width, target_height):
    image = tf.keras.preprocessing.image.load_img(img_file)
    width, height = image.size
    left   = width*0.12
    top    = height*0.12
    right  = width*0.88
    bottom = height*0.88
    image  = image.crop((left, top, right, bottom))
    image  = image.resize((target_width, target_height))
    return image

def loadImgs(cnv_list, img_width, img_height):
    cnv_np = np.zeros((len(cnv_list), img_width, img_height, 3))
    for index, each_cnv in enumerate(cnv_list):
        if index % 100 == 1:
            time_stamp = datetime.datetime.now()
            time_str   = time_stamp.strftime('%Y.%m.%d-%H:%M:%S')
            print("[%s] Processing %d %s..."%(time_str, index, each_cnv))
        try:
            cnv_img = resizeCropImg(each_cnv, img_width, img_height)
            cnv_np[index] = tf.keras.preprocessing.image.img_to_array(cnv_img)
        except:
            cnv_np[index] = None
    time_stamp = datetime.datetime.now()
    time_str   = time_stamp.strftime('%Y.%m.%d-%H:%M:%S')
    print("[%s] Done %d."%(time_str, index))
    return cnv_np

'''
Processing data
'''
def fetch_df_by_wins(data_df, data_img, data_label, min_win, max_win, sort=False):
    data_df = data_df.reset_index(drop=True)
    selected_data_df  = data_df[(data_df['Num_Targets_Wins']>=min_win) & (data_df['Num_Targets_Wins']<=max_win)]
    
    if sort != False:
        selected_data_df = selected_data_df.sort_values(by=['Num_Targets_Wins'])
    selected_index    = selected_data_df.index
    
    selected_data_img = data_img[selected_index]
    selected_label    = data_label[selected_index]
    if selected_index.shape[0] != selected_data_img.shape[0] != len(selected_label):
        print("[Error]")
        pdb.set_trace()
        
    print("There is/are %d CNVs with number of targets/windows between %d and %d."%(len(selected_label),min_win, max_win))
    
    return selected_data_df, selected_data_img, selected_label

def fetch_roc_info_by_num_win(model, cnv_info_df, img, label, min_win, max_win):
    print("Processing num_windows: %d-%d ..."%(min_win, max_win))
    selected_df, selected_img, selected_label = fetch_df_by_wins(cnv_info_df, 
                                                                img, 
                                                                label,
                                                                min_win,
                                                                max_win)
    num_cnv = len(selected_df)
    if min_win != max_win:
        roc_info = "#Targets: " + str(min_win) + "-" + str(max_win) + ", #CNVs: " + str(num_cnv)
    else:
        roc_info = "#Targets: " + str(min_win) + ", #CNVs: " + str(num_cnv)
        
    fpr, tpr, thresholds, auc = pred_roc_data(model, selected_img, selected_label)
    
    return fpr, tpr, auc, roc_info

def recall_m(y_true, y_pred):
    true_positives = keras.backend.sum(keras.backend.round(keras.backend.clip(y_true * y_pred, 0, 1)))
    possible_positives = keras.backend.sum(keras.backend.round(keras.backend.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + keras.backend.epsilon())
    return recall

def precision_m(y_true, y_pred):
    true_positives = keras.backend.sum(keras.backend.round(keras.backend.clip(y_true * y_pred, 0, 1)))
    predicted_positives = keras.backend.sum(keras.backend.round(keras.backend.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + keras.backend.epsilon())
    return precision

def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+keras.backend.epsilon()))
#ref: https://datascience.stackexchange.com/questions/45165/how-to-get-accuracy-f1-precision-and-recall-for-a-keras-model

def pred_roc_data(model, img, label_one_hot):
    pred_keras = model.predict(img).ravel() # ravel(): Flatten the array
    fpr_keras, tpr_keras, thresholds_keras = roc_curve(label_one_hot.ravel(), pred_keras, drop_intermediate=False) 
    auc_keras = metrics.auc(fpr_keras, tpr_keras)
    return fpr_keras, tpr_keras, thresholds_keras, auc_keras

def output_model_metrics(model_name, loss, accuracy, f1_score, precision, recall, output_file):
    model_metric_list = []
    model_metric_list.append([model_name])
    model_metric_list.append(["loss:", loss]) 
    model_metric_list.append(["accuracy:", accuracy])
    model_metric_list.append(["precision:", precision])
    model_metric_list.append(["recall:", recall])
    model_metric_list.append(["f1_score:", f1_score])
    func.output_to_file(model_metric_list, output_file)  

'''
show results
''' 
#from sklearn.metrics import classification_report
def show_confusion_matrix(validations, predictions, labels, output_img_file=None):
    matrix = metrics.confusion_matrix(validations, predictions)
    plt.figure(figsize=(6, 4),dpi=150)
    sns.heatmap(matrix,
                cmap="coolwarm",
                linecolor='white',
                linewidths=1,
                xticklabels=labels,
                yticklabels=labels,
                annot=True,
                fmt="d")
    plt.title("Confusion Matrix", fontsize=16)
    plt.ylabel("True Label", fontsize=16)
    plt.xlabel("Predicted Label", fontsize=16)    
    
    if output_img_file != None:
        plt.savefig(output_img_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("Figure has been output plot to:",output_img_file)
    plt.show()
    plt.close()
    
def confusion_matrix(model, test_img, test_label, nClasses, output_img_file=None):
    print("\n--- Confusion matrix for test data ---\n")
    if nClasses == 4:
        print("4 classes label: 0-True del; 1-False del; 2-False dup; 3-True dup")
        labels = ['DEL', 'Not_DEL', 'Not_DUP', 'DUP']
    if nClasses == 3:
        print("3 classes label: 0-True del; 1-True dip; 2-True dup")
        labels = ['DEL', 'DIP', 'DUP']   
    if nClasses == 2:
        print("2 classes label: 1-True; 0-False")
        labels = ['False', 'True']  

    test_pred = model.predict(test_img)
    if test_label.ndim >1:
        test_label_one_hot = test_label
    else:
        test_label_one_hot = tf.keras.utils.to_categorical(test_label)
    # Take the class with the highest probability from the test predictions
    max_pred_test = np.argmax(test_pred, axis=1)
    max_label_test = np.argmax(test_label_one_hot, axis=1)
    show_confusion_matrix(max_label_test, max_pred_test, labels, output_img_file)

    print("\n--- Classification report for test data ---\n")
    print(classification_report(max_label_test, max_pred_test))
   
def draw_loss_accuracy_curves(history, project_name, output_img_file=None):
    plt.figure(figsize=[8,6])
    plt.plot(history.history['loss'],'r',linewidth=3.0)
    plt.plot(history.history['val_loss'],'b',linewidth=3.0)
    plt.legend(['Training loss', 'Validation Loss'],fontsize=18)
    plt.xlabel('Epochs ',fontsize=16)
    plt.ylabel('Loss',fontsize=16)
    plt.title('Loss Curves',fontsize=16)

    plt.figure(figsize=[8,6])
    plt.plot(history.history['accuracy'],'r',linewidth=3.0)
    plt.plot(history.history['val_accuracy'],'b',linewidth=3.0)
    plt.legend(['Training Accuracy', 'Validation Accuracy'],fontsize=18)
    plt.xlabel('Epochs ',fontsize=16)
    plt.ylabel('Accuracy',fontsize=16)
    plt.ylim(0.7, 1)
    plt.title(project_name+' accuracy curves',fontsize=16)

    if output_img_file != None:
        plt.savefig(output_img_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("Figure has been output plot to:",output_img_file)
    plt.show()
    plt.close()

   
def draw_single_roc_curve(tpr, fpr, auc, output_img_file=None):
    plt.figure(1,dpi=150)
    plt.tick_params(labelsize="x-large")
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr,tpr,label='AUC area = {:.3f})'.format(auc))
    plt.xlabel('False positive rate',fontsize="xx-large")
    plt.ylabel('True positive rate',fontsize="xx-large")
    plt.title('ROC curve',fontsize="xx-large")
    plt.legend(loc='best',fontsize="large")
    
    if output_img_file != None:
        plt.savefig(output_img_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("ROC curve output plot to:",output_img_file)
    plt.show()
    plt.close() 
    
    # Zoom in view of the upper left corner.
    plt.figure(2,dpi=150)
    plt.tick_params(labelsize="x-large")
    plt.xlim(0, 0.3)
    plt.ylim(0.7, 1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr,tpr,label='AUC area = {:.3f})'.format(auc))
    plt.xlabel('False positive rate',fontsize="xx-large")
    plt.ylabel('True positive rate',fontsize="xx-large")
    plt.title('ROC curve (zoomed in at top left)',fontsize="xx-large")
    plt.legend(loc='best',fontsize="large")

    if output_img_file != None:
        path, filename, file_extension = func.extractFilePathNameExtension(output_img_file)
        image_zoom_file = path + '/' + filename + "_zoom" +file_extension
        plt.savefig(image_zoom_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("Zoomed ROC curve output plot to:",image_zoom_file)
    plt.show()
    plt.close()
    

def draw_multiple_roc_curve(tpr_list, fpr_list, auc_list, info_list, title_content='ROC curve', output_image_file=None):   
    label_size = 16 #"x-large"
    
    plt.figure(1,dpi=150)
    plt.tick_params(labelsize= label_size)
    plt.plot([0, 1], [0, 1], 'k--')
        
    for roc_i in range(len(auc_list)):
        plt.plot(fpr_list[roc_i], tpr_list[roc_i], lw=2, alpha=0.3, 
                 label='%s (AUC = %0.3f)' % (info_list[roc_i], auc_list[roc_i]))

    plt.xlabel('False positive rate',fontsize=label_size)
    plt.ylabel('True positive rate', fontsize=label_size)
    plt.title(title_content,fontsize=label_size)
    plt.legend(loc='best',fontsize=label_size/2)
    
    if output_image_file != None:
        plt.savefig(output_image_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("ROC curve output plot to:",output_image_file)
    plt.show()
    plt.close() 
    
    ## draw the zoomed in figure
    plt.figure(2,dpi=150)
    plt.tick_params(labelsize=label_size)
    plt.xlim(0, 0.3)
    plt.ylim(0.7, 1)
    plt.plot([0, 1], [0, 1], 'k--')

    for roc_i in range(len(auc_list)):
        plt.plot(fpr_list[roc_i], tpr_list[roc_i], lw=2, alpha=0.3, 
                 label='%s (AUC = %0.3f)' % (info_list[roc_i], auc_list[roc_i]))

    plt.xlabel('False positive rate',fontsize=label_size)
    plt.ylabel('True positive rate',fontsize=label_size)
    plt.title(title_content ,fontsize=label_size) #'(zoomed in at top left)'
    plt.legend(loc='best',fontsize=label_size/2)
    
    if output_image_file != None:
        path, filename, file_extension = func.extractFilePathNameExtension(output_image_file)
        image_zoom_file = path + '/' + filename + "_zoom" +file_extension
        plt.savefig(image_zoom_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("Zoomed ROC curve output plot to:",image_zoom_file)
    plt.show()
    plt.close()
    
def draw_kfold_roc_curve(tpr_list, tpr_interp_list, fpr_list, auc_list, output_image_file=None, std=False):
    if len(fpr_list) == len(tpr_list) == len(tpr_interp_list) == len(auc_list):
        num_folds = len(auc_list)
        print("Number of folds:",num_folds)
    else:
        print("[Error] please check the folllowing numbers.")
        print("length of tpr_list:",len(tpr_list))
        print("length of tpr_interp_list:",len(tpr_interp_list))
        print("length of fpr_list:",len(fpr_list))
        print("length of auc_list:",len(auc_list))
        
    plt.figure(1,dpi=150)
    plt.tick_params(labelsize=16)
    plt.plot([0, 1], [0, 1], 'k--')
        
    mean_fpr = np.linspace(0, 1, 100)
    for fold_num in range(num_folds):
        plt.plot(fpr_list[fold_num], tpr_list[fold_num], lw=2, alpha=0.3, label='ROC fold %d (AUC = %0.3f)' % (fold_num, auc_list[fold_num]))

    mean_tpr = np.mean(tpr_interp_list, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = sklearn.metrics.auc(mean_fpr, mean_tpr)

    std_auc = np.std(auc_list)
    plt.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    ## plot std
    if std == True:
        std_tpr = np.std(tpr_list, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                        label=r'$\pm$ 1 std. dev.')

    plt.xlabel('False positive rate',fontsize=16)
    plt.ylabel('True positive rate',fontsize=16)
    plt.title('ROC curve',fontsize=16)
    plt.legend(loc='best',fontsize=9)
    if output_image_file != None:
        plt.savefig(output_image_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("ROC curve output plot to:",output_image_file)
    plt.show()
    plt.close() 
    func.showDateTime()
    
    ## draw the zoomed in figure
    plt.figure(2,dpi=150)
    plt.tick_params(labelsize="x-large")
    plt.xlim(0, 0.3)
    plt.ylim(0.7, 1)
    plt.plot([0, 1], [0, 1], 'k--')

    for fold_num in range(num_folds):
        plt.plot(fpr_list[fold_num], tpr_list[fold_num], lw=2, alpha=0.3, label='ROC fold %d (AUC = %0.3f)' % (fold_num, auc_list[fold_num]))

    mean_tpr = np.mean(tpr_interp_list, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = sklearn.metrics.auc(mean_fpr, mean_tpr)
    std_auc = np.std(auc_list)
    plt.plot(mean_fpr, mean_tpr, color='b',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    plt.xlabel('False positive rate',fontsize=16)
    plt.ylabel('True positive rate',fontsize=16)
    plt.title('ROC curve (zoomed in at top left)',fontsize=16)
    plt.legend(loc='best',fontsize=9)
    if output_image_file != None:
        path, filename, file_extension = func.extractFilePathNameExtension(output_image_file)
        image_zoom_file = path + filename + "_zoom" +file_extension
        plt.savefig(image_zoom_file, facecolor='w', edgecolor='w', bbox_inches = 'tight')
        print("Zoomed ROC curve output plot to:",image_zoom_file)
    plt.show()
    plt.close()
   
# ROC curve with color bar to view thresholds
def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates,
    in the correct format for LineCollection:
    an array of the form
    numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    return segments


def colorline(x, y, z=None, axes=None,
              cmap=plt.get_cmap('coolwarm'),
              norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0,
              **kwargs):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if isinstance(z, numbers.Real):
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = matplotlib.collections.LineCollection(
        segments, array=z, cmap=cmap, norm=norm,
        linewidth=linewidth, alpha=alpha, **kwargs
    )

    if axes is None:
        axes = plt.gca()

    axes.add_collection(lc)
    axes.autoscale()

    return lc


def plot_roc(tpr, fpr, thresholds, subplots_kwargs=None,
             label_every=None, label_kwargs=None,
             fpr_label='False Positive Rate',
             tpr_label='True Positive Rate',
             luck_label='Luck',
             title='Receiver operating characteristic',
             **kwargs):

    if subplots_kwargs is None:
        subplots_kwargs = {}

    figure, axes = plt.subplots(1, 1, **subplots_kwargs)

    if 'lw' not in kwargs:
        kwargs['lw'] = 1

    axes.plot(fpr, tpr, **kwargs)

    if label_every is not None:
        if label_kwargs is None:
            label_kwargs = {}

        if 'bbox' not in label_kwargs:
            label_kwargs['bbox'] = dict(
                boxstyle='round,pad=0.5', fc='yellow', alpha=0.5,
            )

        for k in six.moves.range(len(tpr)):
            if k % label_every != 0:
                continue

            threshold = str(numpy.round(thresholds[k], 2))
            x = fpr[k]
            y = tpr[k]
            axes.annotate(threshold, (x, y), **label_kwargs)

    if luck_label is not None:
        axes.plot((0, 1), (0, 1), '--', color='Gray', label=luck_label)

    lc = colorline(fpr, tpr, thresholds, axes=axes)
    figure.colorbar(lc)

    axes.set_xlim([-0.05, 1.05])
    axes.set_ylim([-0.05, 1.05])

    axes.set_xlabel(fpr_label)
    axes.set_ylabel(tpr_label)

    axes.set_title(title)

    axes.legend(loc="lower right")

    return figure, axes

'''
CNN models
'''
def cnn_model(model_name, nClasses, input_shape):
    if model_name == "CNN_model":
        model = Sequential()
        # The first two layers with 32 filters of window size 3x3
        model.add(Conv2D(32, (3, 3), padding='same', activation='relu', input_shape=input_shape))
        model.add(Conv2D(32, (3, 3), activation='relu'))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))

        model.add(Conv2D(64, (3, 3), padding='same', activation='relu'))
        model.add(Conv2D(64, (3, 3), activation='relu'))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))

        model.add(Conv2D(64, (3, 3), padding='same', activation='relu'))
        model.add(Conv2D(64, (3, 3), activation='relu'))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        model.add(Dropout(0.25))

        model.add(Flatten())
        model.add(Dense(512, activation='relu'))
        model.add(Dropout(0.5))
        model.add(Dense(nClasses, activation='softmax'))
        
        #opt = keras.optimizers.Adam(learning_rate=0.01)
        #The learning rate. Defaults to 0.001
        model.compile(optimizer=tf.keras.optimizers.RMSprop(learning_rate=0.0001), loss='categorical_crossentropy',  metrics=['accuracy', f1_m, precision_m, recall_m])
        return model
    
    elif model_name == "MobileNet_v1":
        model = tf.keras.applications.MobileNet(
                    input_shape=None,
                    alpha=1.0,
                    depth_multiplier=1,
                    dropout=0.001,
                    include_top=True,
                    weights="imagenet",
                    input_tensor=None,
                    pooling=None,
                    classes=1000,
                    classifier_activation="softmax"
                )
        ## self added FC layer        
        new_model = keras.models.Sequential()
        new_model.add(model)
        new_model.add(Flatten())
        new_model.add(Dense(512,activation='relu'))
        new_model.add(Dropout(0.5))
        new_model.add(Dense(nClasses, activation='softmax'))

        new_model.compile(loss='categorical_crossentropy', metrics=['accuracy', f1_m, precision_m, recall_m])
        return new_model
    
    elif model_name == 'ResNet50':
        model = tf.keras.applications.ResNet50(
                include_top=True,
                weights="imagenet",
                input_tensor=None,
                input_shape=None,
                pooling=None,
                classes=1000
            )
        ## self added FC layer        
        new_model = Sequential()
        new_model.add(model)
        new_model.add(Flatten())
        new_model.add(Dense(512,activation='relu'))
        new_model.add(Dropout(0.5))
        new_model.add(Dense(nClasses, activation='softmax'))

        new_model.compile(loss='categorical_crossentropy', metrics=['accuracy', f1_m, precision_m, recall_m])
        return new_model
    
    else:
        print("[Error] No this model%s:"%model_name)
        pdb.set_trace()
        
def transfer_learning_model(model_name, nClasses, learning_rate, trainable=False):
    if model_name == "MobileNet_v1":
        base_model = tf.keras.applications.MobileNet(
        weights='imagenet',  # Load weights pre-trained model.
        input_shape=(224, 224, 3),    #input_shape=(224, 224, 3),
        include_top=False)  # Do not include the ImageNet classifier at the top.
        
        print("Model name: %s, nClasses: %d, Learning rate:%f, Trainable: %s"%(model_name, nClasses, learning_rate, trainable))
        base_model.trainable = trainable
        inputs = keras.Input(shape=(224, 224, 3)) #keras.Input(shape=(224, 224, 3))
        # We make sure that the base_model is running in inference mode here,
        # by passing `training=False`. This is important for fine-tuning, as you will
        # learn in a few paragraphs.
        x = base_model(inputs, training=trainable)
        
        # Convert features of shape `base_model.output_shape[1:]` to vectors
        x = keras.layers.GlobalAveragePooling2D()(x)
        # A Dense classifier with a single unit (binary classification)
        outputs = keras.layers.Dense(nClasses,activation='softmax')(x)
        
        model = keras.Model(inputs, outputs)
        #model.compile(loss='categorical_crossentropy', metrics=['accuracy', f1_m, precision_m, recall_m])
        model.compile(optimizer=keras.optimizers.Adam(learning_rate=learning_rate),
                      loss='categorical_crossentropy', metrics=['accuracy', f1_m, precision_m, recall_m])  
        print("Number of base model layers:", len(base_model.layers))
        return model
    
    if model_name == "ResNet50":
        base_model = tf.keras.applications.ResNet50(
        weights='imagenet',  # Load weights pre-trained model.
        input_shape=(224, 224, 3),    #input_shape=(224, 224, 3),
        include_top=False)  # Do not include the ImageNet classifier at the top.
        
        print("Number of base model layers:", len(base_model.layers))
        print("Model name: %s, nClasses: %d, Learning rate:%f, Trainable: %s"%(model_name, nClasses, learning_rate, trainable))
        base_model.trainable = trainable
        inputs = keras.Input(shape=(224, 224, 3)) #keras.Input(shape=(224, 224, 3))
        # We make sure that the base_model is running in inference mode here,
        # by passing `training=False`. This is important for fine-tuning, as you will
        # learn in a few paragraphs.
        x = base_model(inputs, training=trainable)
        
        # Convert features of shape `base_model.output_shape[1:]` to vectors
        x = keras.layers.GlobalAveragePooling2D()(x)
        # A Dense classifier with a single unit (binary classification)
        outputs = keras.layers.Dense(nClasses,activation='softmax')(x)
        
        model = keras.Model(inputs, outputs)
        #model.compile(loss='categorical_crossentropy', metrics=['accuracy', f1_m, precision_m, recall_m])
        model.compile(optimizer=keras.optimizers.Adam(learning_rate=learning_rate),
                      loss='categorical_crossentropy', metrics=['accuracy', f1_m, precision_m, recall_m]) 
        return model
    
    else:
        print("[Error] No this model%s:"%model_name)
        pdb.set_trace()
