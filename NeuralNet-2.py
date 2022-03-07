#####################################################################################################################
#   Assignment 2: Neural Network Analysis
#   This is a starter code in Python 3.6 for a neural network.
#   You need to have numpy and pandas installed before running this code.
#   You need to complete all TODO marked sections
#   You are free to modify this code in any way you want, but need to mention it
#       in the README file.
#
#####################################################################################################################

#https://archive.ics.uci.edu/ml/datasets/Haberman%27s+Survival
#https://archive.ics.uci.edu/ml/datasets/Mammographic+Mass
import numpy as np
import pandas as pd
import sys
import os
import sklearn
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report


class NeuralNet:
    def __init__(self, dataFile, header=True):
        self.raw_input = pd.read_csv(dataFile)
        

    

    # TODO: Write code for pre-processing the dataset, which would include
    # standardization, normalization,
    #   categorical to numerical, etc
    def preprocess(self):
        self.processed_data = self.raw_input
        
        
        self.processed_data.columns = ['Age','Year','Axillary_node','survival']
        
        #describe data and plot variables 
        print(type(self.processed_data))
        print(self.processed_data.head())
        print(self.processed_data.describe())
        
        
        ##Plot the distribution of variable 
        self.processed_data.hist(figsize=(20,15))
        plt.tight_layout()
        plt.savefig("Histogram.jpg", dpi=150)
        plt.show()
        
        #pair plot
        sns.pairplot(self.processed_data, hue="survival")
        plt.tight_layout()
        plt.savefig("pairplot.jpg", dpi=150)
        plt.show()
        
        #check and remove the missing values
        self.processed_data.isnull().sum() 
        
        # handle the duplicat rows 
        print(self.processed_data.duplicated())
        self.processed_data.drop_duplicates() ## keep the first occurance
        
        # Handle mising values 
        self.processed_data.dropna(axis=0,inplace=True) # Due to the nature of this data,
                                                            #if there is a missing value, the row should be dropped.
                                                            
        
        return 0

    # TODO: Train and evaluate models for all combinations of parameters
    # specified. We would like to obtain following outputs:
    #   1. Training Accuracy and Error (Loss) for every model
    #   2. Test Accuracy and Error (Loss) for every model
    #   3. History Curve (Plot of Accuracy against training steps) for all
    #       the models in a single plot. The plot should be color coded i.e.
    #       different color for each model

    def train_evaluate(self):
        
        ncols = len(self.processed_data.columns)
        nrows = len(self.processed_data.index)
        
        X = self.processed_data.iloc[:, 0:(ncols - 1)]
        y = self.processed_data.iloc[:, (ncols-1)]
        X_train, X_test, y_train, y_test = train_test_split(
            X, y,test_size=0.2, random_state=42)
        
      #TO DO  #normalization/standarization 
        sc = StandardScaler()
        
        #fit the training data
        sc.fit(X_train)
        
        #scaling
        X_train=sc.transform(X_train)
        X_test=sc.transform(X_test)  
        
        # Create the neural network and be sure to keep track of the performance
        # Below are the hyperparameters that you need to use for model
        #   evaluation
        activations = ['logistic', 'tanh', 'relu']
        learning_rate = [0.01, 0.1]
        max_iterations = [100, 200] # also known as epochs
        #num_hidden_layers = [2, 3]
        hidden_layers1 = [(2,2),(3,3)]
        hidden_layers2 = [(2,2,2),(3,3,3)]
        res_li_hl1 = {}
        res_li_hl2 = {}
        for act in activations:
            for lr in learning_rate:
                for itr in max_iterations:
                    for hl in hidden_layers1:
                        cl_h2 = MLPClassifier(hidden_layer_sizes=hl,max_iter=itr,activation=act,
                                learning_rate='constant',
                                learning_rate_init=lr)
                        r1 = cl_h2.fit(X_train,y_train)
                        y_test_pred = cl_h2.predict(X_test)
                        y_train_pred = cl_h2.predict(X_train)
                        res_li_hl1[r1] = [y_test_pred,y_train_pred]
                    
                    for hl2 in hidden_layers2:
                        cl_h3 = MLPClassifier(hidden_layer_sizes=hl2,max_iter=itr,activation=act,
                                learning_rate='constant',
                                learning_rate_init=lr)
                        
                        r2 = cl_h3.fit(X_train,y_train)
                        y_test_pred_h3 = cl_h3.predict(X_test)
                        y_train_pred_h3 = cl_h3.predict(X_train)
                        res_li_hl2[r2] = [y_test_pred_h3,y_train_pred_h3]
                    
                     
                    
        return res_li_hl1,res_li_hl2,y_train,y_test,X_train
   
    #TOUTPUT SUMMARIZATION  
   # 1. Training Accuracy and Error (Loss) for every model
   # 2. Test Accuracy and Error (Loss) for every model
   # 3. History Curve (Plot of loss /Accuracy against training steps) for all
   #     the models in a single plot. The plot should be color coded i.e.
   #     different color for each model
    def summarize_results(self):
        li_dict,li_dicth3,y_train_lab,y_test_lab,X_train = self.train_evaluate()
        
        #summarizing hiddden layer size 2 [(2,2),(3,3)]
        labels_h2 = ['logistic_0.01_h2_2_it100','logistic_0.01_h3_3_it100','logistic_0.01_h2_2_it200',
                  'logistic_0.01_h3_3_it200','logistic_0.1_h2_3_it100','logistic_0.1_h3_3_it100',
                  'logistic_0.1_h2_2_it200','logistic_0.1_h3_3_it200',
                  'tanh_0.01_h2_2_it100','tanh_0.01_h3_3_it100','tanh_0.01_h2_2_it200',
                  'tanh_0.01_h3_3_it200','tanh_0.1_h2_2_it100','tanh_0.1_h3_3_it100',
                  'tanh_0.1_h2_2_it200','tanh_0.1_h3_3_it200',
                  'relu_0.01_h2_2_it100','relu_0.01_h3_3_it100','relu_0.01_h2_2_it200',
                  'relu_0.01_h3_3_it200','relu_0.1_h2_2_it100','relu_0.1_h3_3_it100',
                  'relu_0.1_h2_2_it200','relu_0.1_h3_3_it200']
         
        labels_h3 = ['logistic_0.01_h2_2_2_it100','logistic_0.01_h3_3_3_it100','logistic_0.01_h2_2_2_it200',
                  'logistic_0.01_h3_3_3_it200','logistic_0.1_h2_3_3_it100','logistic_0.1_h3_3_3_it100',
                  'logistic_0.1_h2_2_2_it200','logistic_0.1_h3_3_3_it200',
                  'tanh_0.01_h2_2_2_it100','tanh_0.01_h3_3_3_it100','tanh_0.01_h2_2_2_it200',
                  'tanh_0.01_h3_3_3_it200','tanh_0.1_h2_2_2_it100','tanh_0.1_h3_3_3_it100',
                  'tanh_0.1_h2_2_2_it200','tanh_0.1_h3_3_3_it200',
                  'relu_0.01_h2_2_2_it100','relu_0.01_h3_3_3_it100','relu_0.01_h2_2_2_it200',
                  'relu_0.01_h3_3_3_it200','relu_0.1_h2_2_2_it100','relu_0.1_h3_3_3_it100',
                  'relu_0.1_h2_2_2_it200','relu_0.1_h3_3_3_it200']
                 
                  
        #score_train=[]
        for ky,vl in li_dict.items():
            test_accuracy = accuracy_score(y_test_lab,vl[0])
            train_accuracy = accuracy_score(y_train_lab,vl[1])
            with open('summary_table.hidden_layers2.txt', 'a') as f:
                
            
                print(('model: %s\tTrain accuracy: %.3f\tTest accuracy:%.3f\tMSE Train: %.3f\tMSE Test: %.3f\t'
                  % (ky,train_accuracy,test_accuracy,mean_squared_error(y_train_lab,vl[1]),
                     mean_squared_error(y_test_lab,vl[0]))),file=f)
            
            
            print(classification_report(y_test_lab,vl[0]))
            print("\n")
            print("\n")
        
            
            plt.plot(ky.loss_curve_)
            
            plt.xlabel('Iterations')
            plt.ylabel('cost')
            plt.title("Loss curve (Hidden layer 2)")
            plt.legend(labels=labels_h2,fontsize=4)
            plt.tight_layout()
            plt.savefig(("Loss_curve_hidden_layer_2)"+".jpg"),dpi=300)
            #plt.show()
          
        #summarizing hidden layer 3 , [(2,2,2),(3,3,3)]
        
        for ky_h3,vl in li_dicth3.items():
            test_accuracy_h3 = accuracy_score(y_test_lab,vl[0])
            train_accuracy_h3 = accuracy_score(y_train_lab,vl[1])
            with open('summary_table.hidden_layes_three.txt', 'a') as f:
                
            
                print(('model: %s\tTrain accuracy: %.3f\tTest accuracy:%.3f\tMSE Train: %.3f\tMSE Test: %.3f\t'
                  % (ky_h3,train_accuracy_h3,test_accuracy_h3,mean_squared_error(y_train_lab,vl[1]),
                     mean_squared_error(y_test_lab,vl[0]))),file=f)
            
            
            print(classification_report(y_test_lab,vl[0]))
            print("\n")
            print("\n")
        
            
            plt.plot(ky_h3.loss_curve_)
            
            plt.xlabel('Iterations')
            plt.ylabel('cost')
            plt.title("Loss curve (Hidden layer 3)")
            plt.legend(labels=labels_h3,fontsize=4)
            plt.tight_layout()
            plt.savefig(("Loss_curve_hidden_layer_3"+".jpg"),dpi=300)
            #plt.show()
        
    
        return 0




if __name__ == "__main__":
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/haberman/haberman.data"
    #url = "/Users/ashwani/Documents/ML/Neural_network/haberman.data"
    neural_network = NeuralNet(url) # put in path to your file
    neural_network.preprocess()
    neural_network.train_evaluate()
    neural_network.summarize_results()
