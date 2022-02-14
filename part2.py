#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ashwani
"""

#Import packages

import sys
import os
import numpy as np
import pandas as pd
import sklearn
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.linear_model import SGDRegressor
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score




class LinearRegGD_Implementation():
    def __init__(self,file_path):
        self.file_path = file_path
        
        
    def preprocessing_data(self):
        
        #Read the sensor data
        df = pd.read_csv(self.file_path)
        
        #Shows the distribution of data
        print(df.head())
        print(df.info())
        print(df.describe())
        
        #Check for missing or null values. There was no missing values, hence nothing no columns were removed.
        print (df.isnull().sum())
        
        #Histogram of each variable
        df.hist(bins=50,figsize=(20,15))
        plt.title("Histogram of each independent variable (Distribution of data)")
        plt.savefig("Histogram.jpg", dpi=150)
        plt.show()
        
        
        #Peasron correlation heatmap of features using z score. Here z score is
        #considered just for scaling and visulaiztion
        stsc = StandardScaler()
        df_zscore = stsc.fit_transform(df) # performs standarization by centering and scaling
        #print(type(df_zscore))
        cols = ['anchor_ratio','trans_range','node_density','iterations']
        cor = np.corrcoef(df[cols].values.T)
        sns.set(font_scale=1.5)
        
        
        hm = sns.heatmap(cor,cbar=True,
                         annot=True,
                         square=True,
                         fmt='.2f',
                         annot_kws={'size':15},
                         yticklabels=cols,
                         xticklabels=cols)
        plt.savefig("Correaltion_heatmap.jpg",dpi=100)
        plt.show()
        
        #split the training and testing data.
        X,y = df.iloc[:,:4].values, df.iloc[:,4].values
        x_train,x_test,y_train,y_test =train_test_split(X,y,test_size=0.2,random_state=42)
        
        
        self.predict(x_train,x_test,y_train,y_test)
        
        
        
    def predict(self,X_train,X_test,Y_train,Y_test):
          
        #Print to see the dimension of split data arrays.
        print(X_test.shape)
        print(Y_test.shape)
        print(X_train.shape)
        print(Y_train.shape)
        
        
        #Aplha here is constant that multiplies the regularization term.Also used to compute the learning rate
        #alpha_rate = [0.1,0.01,0.00001,1.28,1.85,3.0] 
        
        rate ='invscaling'
        #eta = eta0 / pow(t, power_t)
        
        #redirect the log file to a text file written in working dir
        sys.stdout = open('logs.%s.txt' % rate ,'w') 
        
        #Transform and fit using pipeline. Test data is only transformed.
        pipe_lr = Pipeline([('scl',StandardScaler()),
                           ('clf',SGDRegressor(max_iter=1000,penalty=None,
                                               learning_rate='invscaling',
                                               shuffle=False,
                                               verbose=1))])
        
        #Fit
        pipe_lr.fit(X_train,Y_train)
        y_train_pred = pipe_lr.predict(X_train)
        y_test_pred = pipe_lr.predict(X_test)
        
        ##Plot the results
        plt.scatter(y_test_pred, y_test_pred - Y_test, c='red',marker='o',
                    label="Test data")
        plt.scatter(y_train_pred, y_train_pred - Y_train, c='blue',marker='s',
                    label="Training data")
        plt.legend()
        plt.xlabel('Predicted value')
        plt.ylabel('Residuals')
        plt.hlines(y=0,xmin=0.25,xmax=1.50,lw=2,color='red')
        plt.title("Residual vs Predicted for: %s "  % (rate))
        plt.savefig(("Residulas_vs_predicted"+str(rate)+".jpg"),dpi=300)
        plt.show()
        
        
        # Evaluating the performance of linear regression model
        print('Learning Rate, %s, MSE Test:, %.3f, MSE Train: %.3f:' % (rate, mean_squared_error(Y_test,y_test_pred),
                                                    mean_squared_error(Y_train,y_train_pred)))      
        
        print('Learning Rate, %s, R^2 Test: %.3f, R^2 Train %.3f' %(rate, r2_score(Y_test,y_test_pred),
                                                 r2_score(Y_train,y_train_pred)))
        
        
        #END
            
        
        
if __name__ == '__main__':
    url = "https://raw.githubusercontent.com/Shashank425/ML/main/mcs_ds_edited_iter_shuffled.csv"
    obj1 = LinearRegGD_Implementation(url)
    obj1.preprocessing_data()
    
