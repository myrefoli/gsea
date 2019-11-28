import sys
import math
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt

class KNN():
    """
    KNN classifies samples as healthy or patient using K-Nearest Neighbors
    """
    def __init__(self):
        self.exp_file = pd.DataFrame([])
        self.samp_file = pd.DataFrame([])

    def load_data(self, exp_file, samp_file):
        """
        Takes the paths to the expression and sample files, reads them in, and stores them within the KNN
        """
        self.exp_file = pd.read_csv(exp_file, sep = '\t')
        self.exp_file = self.exp_file.drop(columns = 'SYMBOL')
        self.samp_file = pd.read_csv(samp_file, sep = '\t', header = None)
        self.samp_file = self.samp_file.values.tolist()
        self.samp_file = dict(self.samp_file)

    def get_assignments(self, k, fn):
        """
        Given k neighbors and fn fraction for positive classification, returns class assignments for all samples

        Inputs:
            k is an int representing the given number of neighbors
            fn is a float representing the given fraction needed for a positive classification
        Returns:
            a list of integer 0s and 1s representing class assignments for all samples
        """
        unsorted_assignments = []
        samp_names = []

        #LOOCV for all samples
        for samp_name, samp_val in self.exp_file.items():
            #creating dict of samples
            samp_names.append(samp_name)
            
            temp_exp_file = self.exp_file
            temp_exp_file = temp_exp_file.drop(columns = samp_name)
            num_positives = 0
            distances = []

            #computing euclidean distances
            total_distance = 0
            for comp_samp_name, comp_samp_val in temp_exp_file.items():
                all_distances = samp_val - comp_samp_val
                
                for distance in all_distances:
                    total_distance += distance**2
                
                total_distance = math.sqrt(total_distance)
                distances.append((total_distance, comp_samp_name))

            #sorting for k closest samples
            distances.sort()

            #assigning as patient or healthy
            neighbors = 0
            for i in range(k):
                curr_samp = distances[i]
                curr_samp_name = curr_samp[1]
                samp_class = self.samp_file[curr_samp_name]

                if samp_class == 1:
                    neighbors += 1
            
            if neighbors > (fn * k):
                unsorted_assignments.append(1)

            else:
                unsorted_assignments.append(0)           

        #sorting to match order in samples file
        exp_test = dict(zip(samp_names, unsorted_assignments))
        
        class_assignments = []
        for samp_name in self.samp_file:
            class_assignments.append(exp_test[samp_name])

        return class_assignments

    def calc_metrics(self, k, fn):
        """
        Given k neighbors and fn fraction for positive classification, returns sensitivity and specificity of a KNN classifier

        Inputs: 
            k is an int representing the given number of neighbors
            fn is a float representing the given fraction needed for a positive classification
        Returns:
            a list of floats [sensitivity, specificity] of a KNN classifier using k and fn
        """
        assignments = self.get_assignments(k, fn)
        labels = []
        metrics = []
        id_pos = 0
        all_pos = 0
        id_neg = 0
        all_neg = 0
        
        #getting actual labels of samples
        for samp_name in self.samp_file:
            labels.append(self.samp_file[samp_name])
        
        #computing all positives and negatives
        for label in labels:
            if label == 0:
                all_neg += 1
            else:
                all_pos += 1

        #computing identified positives and negatives
        for i in range(len(assignments)):
            if assignments[i] == 0 and assignments[i] == labels[i]:
                id_neg += 1

            if assignments[i] == 1 and assignments[i] == labels[i]:
                id_pos += 1

        metrics.append(id_pos / all_pos)
        metrics.append(id_neg / all_neg)
        print(metrics)
        return metrics

    """
    def get_roc_curve(self, k, fn):
        labels = []
        classification = self.get_assignments(k, fn)
        
        for gene in self.samp_file:
            labels.append(self.samp_file[gene])

        print(labels)
        print(classification)
    """
        

def main():

    #check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an expression file and a sample file as args.")
        return
        
    #input variables
    exp_file = sys.argv[1]
    samp_file = sys.argv[2]

    knn = KNN()
    knn.load_data(exp_file, samp_file)
    #knn.get_assignments(6, 0.75)
    #knn.calc_metrics(5, 0.75)
    #sense = []
    #spec = []
    labels = []
    
    for sample in knn.samp_file:
        labels.append(knn.samp_file[sample])

    for k in [1, 2, 3, 4, 5, 6]:
        classification = knn.get_assignments(k, 0.5)
        print(accuracy_score(labels, classification))

if __name__=="__main__":
    main()