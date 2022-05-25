from __future__ import print_function
import os
import sys
import pdb
import copy
import sklearn
import tensorflow as tf
import pandas as pd
import PIL
import numpy as np
import keras
from keras.models import load_model
import function as func
import function_dl as func_dl

def cnn_prediction(cnv_file, model_file, use_gpu, output_file):
    # GPU or CPU selection
    if use_gpu == False or use_gpu == 'False':
        print("Using CPU ...")

    else:
        cuda_available = tf.test.is_built_with_cuda()
        gpu_availabel  = tf.config.list_physical_devices('GPU')
        os.environ['CUDA_VISIBLE_DEVICES'] = "-1"
        if cuda_available and gpu_availabel:
            print("Using GPU ...")
            physical_devices = tf.config.experimental.list_physical_devices('GPU') 
        else:
            print("GPU is not available. Try CPU ...")
            physical_devices = tf.config.experimental.list_physical_devices('CPU')

    # Initial variables
    img_width, img_height = 224, 224

    ## Load pre-calculated model
    custom_objects = {"f1_m":func_dl.f1_m, "precision_m":func_dl.precision_m, "recall_m":func_dl.recall_m}

    model_name = 'MobileNet_v1'
    func.showDateTime('\t')
    print("Loading %s ... from %s"%(model_name, model_file))
    try:
        MobileNet_model = keras.models.load_model(model_file, custom_objects=custom_objects)
        print("Model Loaded. ")
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise    

    # Loading CNV_info and images. 
    func.showDateTime('\t')
    print("Loading CNV info and images ...")
    cnv_info_df = pd.read_csv(cnv_file)
    entire_cnv_images_path_list = cnv_info_df['img_path']
    CNV_TYPE_list = func.global_variables()['CNV_TYPE']
    CNV_TYPE      = func.fetch_colName(CNV_TYPE_list,cnv_info_df.columns)[1]
    img_np        = func_dl.loadImgs(entire_cnv_images_path_list, img_width, img_height)

    # ## Normalization
    # Find the shape of input images and create the variable input_shape
    nRows,nCols,nDims = img_np.shape[1:]
    input_shape       = (nRows, nCols, nDims)
    print("The shape of input tensor:",input_shape)

    # Change to float datatype
    img_np = img_np.astype('float32')

    # Scale the data to lie between 0 to 1
    img_np /= 255

    # Find the unique numbers from the train labels
    nClasses = 3
    print('Total number classes: ', nClasses)

    # Prediction 
    print("Predict the type of copy number images by CNV-espresso ...")
    img_pred = MobileNet_model.predict(img_np)
    pred_output_df = copy.deepcopy(cnv_info_df)
    pred_output_df.insert(pred_output_df.shape[1], 'Prob_DEL', "-")
    pred_output_df.insert(pred_output_df.shape[1], 'Prob_Artifacts', "-")
    pred_output_df.insert(pred_output_df.shape[1], 'Prob_DUP', "-")
    pred_output_df.insert(pred_output_df.shape[1], 'Prob', "-")
    pred_output_df.insert(pred_output_df.shape[1], 'Prediction', "-")
    pred_output_df.insert(pred_output_df.shape[1], 'Status', "-")

    num, correct_count = 0, 0
    for i in range(len(img_pred)):
        num += 1
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prob_DEL')] = img_pred[i][0]
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prob_Artifacts')] = img_pred[i][1]
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prob_DUP')] = img_pred[i][2]
        pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prob')] = np.max(img_pred[i])
            
        if  np.any(np.isnan(img_pred[i])) == True:
            pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] = "NaN"
        else:
            if np.argmax(img_pred[i]) == 0:
                pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] = "DEL"
            elif np.argmax(img_pred[i]) == 1:
                pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] = "Artifacts"
            elif np.argmax(img_pred[i]) == 2:
                pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] = "DUP"
            else:
                pdb.set_trace()
            
        pred_output_df_type_col = func.fetch_colName(CNV_TYPE_list,pred_output_df.columns)[0]
        if pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] == "NaN":
            pred_output_df.iloc[i,pred_output_df.columns.get_loc('Status')] = "Error"
        elif pred_output_df.iloc[i,pred_output_df.columns.get_loc('Prediction')] == pred_output_df.iloc[i, pred_output_df_type_col].upper():
            pred_output_df.iloc[i,pred_output_df.columns.get_loc('Status')] = "Positive"
        else:
            pred_output_df.iloc[i,pred_output_df.columns.get_loc('Status')] = "Negative"

    ## output to file
    pred_output_df.to_csv(output_file,index=False)
    func.showDateTime()
    print("Done.")
