# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 18:39:07 2020

@author: Stavros
"""

import numpy as np
import numpy
import sys
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn import datasets
from sklearn.model_selection import train_test_split
import random
from scipy.stats import norm

def plot_decision_boundary(X, y ,Classifier):
    '''
	Plot the decision boundary of a trained classifier for the given data
	
	Parameters
	----------
	X : 2-d matrix
	    samples observation
	y : array_like
	    class labels for each sample. It must have the same row length as X
	Classifier : classifier object
	    classifier object previously trained (i.e. fitted to the training data
        using the fit(X,y) method).
	Returns
	-------
	None
	'''
    if type(Classifier) == type(NearestCentroid()):
        plot_decision_boundary(X,y,
                    KNeighborsClassifier(n_neighbors=1).fit(Classifier.centroids_,np.unique(y)))
        return
    classes = np.unique(y)
    plt.figure()
    plt.plot(X[y==classes[0],0],X[y==classes[0],1],'o',color='blue')
    plt.plot(X[y==classes[1],0],X[y==classes[1],1],'*',color='red')

    nx, ny = 200, 100
    x_min, x_max = plt.xlim()
    y_min, y_max = plt.ylim()
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, nx),np.linspace(y_min, y_max, ny))
    Z = Classifier.predict_proba(np.c_[xx.ravel(), yy.ravel()])
    Z = Z[:, 1].reshape(xx.shape)
    plt.contour(xx, yy, Z, [0.5], linewidths=2., colors='k')
    plt.show()
    
def learning_curve():
    '''
	Plot the learning curve of the k-Nearest neighbor classifier using two 
    gaussian shaped data while changing the training set size.
    
	Parameters
	----------
	None
    
	Returns
	-------
	None
	'''
    Gaussian_shaped = datasets.make_blobs(n_samples=1000,n_features=2,
                                      centers=[[1,1],[2,2]],cluster_std=[1,2])
    data = Gaussian_shaped[0]
    labels = Gaussian_shaped[1]
    Classifier = KNeighborsClassifier(n_neighbors=1)
    X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.5)
    Sampling_prop = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    error_mean = np.zeros(len(Sampling_prop))
    error_std = np.zeros(len(Sampling_prop))
    
    for j in range(0,len(Sampling_prop)):
        error = np.zeros(1000)
        for i in range(0,1000):
            idx = random.sample(range(0,len(y_train)),int(Sampling_prop[j]*len(y_train)))
            Classifier.fit(X_train[idx,:],y_train[idx])
            y_pred = Classifier.predict(X_test)
            error[i] = sum(y_test != y_pred)/len(y_pred) * 100
        error_mean[j] = np.mean(error)
        error_std[j] = np.std(error)
    
    plt.figure()
    plt.errorbar(Sampling_prop*100/2,error_mean,yerr=error_std)
    plt.title('Learning Curve', fontsize=25)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.xlabel('Training set size %',fontsize=20)
    plt.ylabel('Classification error %',fontsize=20)
    plt.show()
    

def ztest_1samp(a, popmean, sigma, alternative='two-sided', axis=0): 
	'''
	Calculates the Z-test for the mean of ONE group of scores.
	
	Parameters
	----------
	a : array_like
	    sample observation
	popmean : float or array_like
	    expected value in null hypothesis. If array_like, then it must have the
	    same shape as `a` excluding the axis dimension
	sigma : float or array_like
	    the standard deviation of the data. If array_like, then it must have the
	    same shape as `a` excluding the axis dimension
	alternative: less', 'two-sided', or 'greater'
            Whether to get the p-value for the one-sided hypothesis ('less'
            or 'greater') or for the two-sided hypothesis ('two-sided').
	axis : int or None, optional
	    Axis along which to compute test. If None, compute over the whole
	    array `a`.
	Returns
	-------
	statistic : float or array
	    z-statistic
	pvalue : float or array
	    p-value
	'''
	m = np.mean(a, axis=axis)
	N = a.shape[0]
	s = sigma / np.sqrt(N)

	z = (m - popmean) / s

	if alternative == 'two-sided':
		p = 2.0 * norm.sf(np.abs(z)) 


	elif alternative == 'greater':
		p = norm.sf(np.abs(z)) 

	elif alternative == 'less':
		p = norm.cdf(np.abs(z))

	else:
		raise ValueError ('Alternative has to be \'two-sided\', \'greater\' or \'less\'')

	return (z, p)


def ztest_rel(a, b, sigma, alternative='two-sided', axis=0): 
	'''
	Calculates the Z-test for the mean of two PAIRED groups of scores.
	
	Parameters
	----------
	a : array_like
	    group 1 observations
	b : array_like
	    group 2 observations
	sigma : float or array_like
	    the standard deviation of the difference. If array_like, then it must have the
	    same shape as `a` excluding the axis dimension
	alternative: less', 'two-sided', or 'greater'
            Whether to get the p-value for the one-sided hypothesis ('less'
            or 'greater') or for the two-sided hypothesis ('two-sided').
	axis : int or None, optional
	    Axis along which to compute test. If None, compute over the whole
	    array `a`.
	Returns
	-------
	statistic : float or array
	    z-statistic
	pvalue : float or array
	    p-value
	'''
	if a.shape != b.shape:
		raise ValueError ('a and b must be of the same size')

	m = np.mean(a - b, axis=axis)
	N = a.shape[0]
	s = sigma / np.sqrt(N)

	z = m / s

	if alternative == 'two-sided':
		p = 2.0 * norm.sf(np.abs(z)) 


	elif alternative == 'greater':
		p = norm.sf(np.abs(z)) 

	elif alternative == 'less':
		p = norm.cdf(np.abs(z))

	else:
		raise ValueError ('Alternative has to be \'two-sided\', \'greater\' or \'less\'')

	return (z, p)


def ztest_ind(a, b, sigma_a, sigma_b, alternative='two-sided', axis=0): 
	'''
	Calculates the Z-test for the mean of two INDEPENDENT groups of scores.
	
	Parameters
	----------
	a : array_like
	    group 1 observations
	b : array_like
	    group 2 observations
	sigma_a : float or array_like
	    the standard deviation of the first group. If array_like, then it must have the
	    same shape as `a` excluding the axis dimension
	sigma_b : float or array_like
	    the standard deviation of the second group. If array_like, then it must have the
	    same shape as `a` excluding the axis dimension
	alternative: less', 'two-sided', or 'greater'
            Whether to get the p-value for the one-sided hypothesis ('less'
            or 'greater') or for the two-sided hypothesis ('two-sided').
	axis : int or None, optional
	    Axis along which to compute test. If None, compute over the whole
	    array `a`.
	Returns
	-------
	statistic : float or array
	    z-statistic
	pvalue : float or array
	    p-value
	'''
	if a.shape != b.shape:
		raise ValueError ('a and b must be of the same size')

	m = np.mean(a - b, axis=axis)
	N1 = a.shape[0]
	N2 = b.shape[0]
	s = np.sqrt(((sigma_a ** 2) /  N1) + ((sigma_b ** 2) /  N2))

	z = m / s

	if alternative == 'two-sided':
		p = 2.0 * norm.sf(np.abs(z)) 

	elif alternative == 'greater':
		p = norm.sf(np.abs(z)) 

	elif alternative == 'less':
		p = norm.cdf(np.abs(z))

	else:
		raise ValueError ('Alternative has to be \'two-sided\', \'greater\' or \'less\'')

	return (z, p)


#Takes as input a distance matrix of pairwise distances and a list of point labels and returns a dictionary with two point labels delimited by comma as key and the distance between those points as value. 
def annotate_distances(dist_mat,label_list):
    points = {}
    index = 0
    for i in range(len(label_list)):
        for j in range(i+1,len(label_list)):
                points[label_list[i]+','+label_list[j]] = dist_mat[index]
                index +=1
    return points


#Takes as input a nx2 numpy array and optional a list of labels of the size n. The function returns NULL and creates a scatter plot with un-/labelled points. 
def plot_scatter(input_dat,label_list=[]):
    for i in range(len(input_dat)):
        x = input_dat[i][0]
        y = input_dat[i][1]
        plt.plot(x, y, 'bo')
        if label_list!=[]:
            plt.text(x * (1 + 0.01), y * (1 + 0.01) , label_list[i], fontsize=12)
    
    plt.xlim((input_dat.min()-0.5, input_dat.max()+0.5))
    plt.ylim((input_dat.min()-0.5, input_dat.max()+0.5))
    plt.show()



def kmeans_cluster_vis(input_data=None,labels=None,cluster_centers=None,dimensions=(0,1)):
    '''
    Plots the clusters and cluster centers as a 2-dimensional scatterplot.
    
    Parameters
    ----------
    input_data : array_like
        The original data array to which clusters were fitted. The array must have objects/samples per row and dimensions per column.
    labels : array_like
        A list or 1-dimensional array with the cluster label for each data point in the input array.
    cluster_centers : array_like
        An array containing the coordinates of the cluster centers. Each row contains 1 cluster center. Must have the same dimensions as the data set.
    dimensions : tuple of ints
        The tuple of integers is indicating which dimensions of the data set will be plotted. Counting dimensions starts with 0. The first integer indicates the data dimension plotted on the x-axis and the second one plotted is plotted on the y-axis.
        default: 0,1
    Returns
    -------
    None
    In terminal, plotted 2-dimensional scatter plot with cluster means and data points colored with respect to the cluster they belong to.
    '''

    #############################
    #checking of input parameters
    #############################
    
    #checking if any parameter is missing
    if input_data is None:
        print("The original data set that was used to fit the KMeans model is missing.")
        return None
    if labels is None:
        print("The cluster labels are missing")
        return None 
    if cluster_centers is None:
        print("The cluster centers are missing")
        return None
    
    #converting parameters to numpy arrays if not already the case
    if type(input_data)!=numpy.ndarray:
        input_data = numpy.array(input_data)
    if type(labels)!=numpy.ndarray:
        labels = numpy.array(labels)
    if type(cluster_centers)!=numpy.ndarray:
        cluster_centers = numpy.array(cluster_centers)
    
    #checking if labels has as many unique values as there are cluster centers
    if len(numpy.unique(labels))!=len(cluster_centers):
        print("The number of provided unique cluster labels and the number of clusters differ even though both should be of size k. \n Number unique cluster labels: %s\n Number of clusters: %s"%(len(numpy.unique(labels)),len(cluster_centers)),file=sys.stderr)
        return None
    
    #checking if labels is 1-dimensional.
    if len(labels.shape)!=1:
       print("The provided cluster labels are a multi dimensional array even though it has to be a list or 1-dimensional array. \n Number of dimensions in cluster labels: %s"%(len(labels.shape)),file=sys.stderr)
       return None
   
    #checking if there are as many labels as samples in the data set.
    if len(input_data)!=len(labels):
        print("The data set and the cluster labels are of different size even though there has to be one label for each data point of the data set. \n Number of samples in data set: %s\nNumber of labels: %s"%(len(input_data),len(labels)),sys.stderr)
        return None
    
    #checking if the cluster_centers and the data set have the same number of dimensions.
    if cluster_centers.shape[1]!=input_data.shape[1]:
        print("The cluster centers and the input data have different numbers of dimensions which must not be. \n Number of cluster center dimensions: %s\n Number of data set dimensions: %s"%(cluster_centers.shape[1],input_data.shape[1]),file=sys.stderr)
        return None
    
    #checking if the data set and cluster_centers have at least as many dimensions as the ones choosen to be plotted.
    if (numpy.max(dimensions)>cluster_centers.shape[1]-1):
        print("The requested dimension for plotting is larger than what the cluster centers or the data set provide. Choose smaller dimensions for plotting.\n Requested dimensions for plotting: %s\n Number of dimensions of the Data set (0-based): %s"%(dimensions,cluster_centers.shape[1]-1),file=sys.stderr)
        return None

    #Checking if the selected number of clusters is within the predefined margins    
    colour = ['b','g','r','c','m','y','w']
    marker = ['o','D','v','^','s','p','*','H','+']
    if len(cluster_centers) > (len(colour)*len(marker)):
        print("You chose a too large k. Maximum k is %s"%(len(colour)*len(marker)),file=sys.stderr)
        return None
    else:
        cluster_colour = []
        for i in marker:
            for j in colour:
                cluster_colour.append(i+j)
    
    #############################
    #Plotting
    #############################

    #Plotting data points
    for i in range(len(input_data)):
        x = input_data[i][dimensions[0]]
        y = input_data[i][dimensions[1]]
        plt.plot(x, y, cluster_colour[labels[i]])
    
    #Plotting cluster centers
    for i in range(len(cluster_centers)):
        x = cluster_centers[i][dimensions[0]]
        y = cluster_centers[i][dimensions[1]]
        plt.plot(x, y, 'kx')
        plt.text(x * (1 + 0.01), y * (1 + 0.01) , "Cluster %s"%str(i), fontsize=12)
        
    plt.xlim((input_data[:,dimensions[0]].min()-0.5, input_data[:,dimensions[0]].max()+0.5))
    plt.ylim((input_data[:,dimensions[1]].min()-0.5, input_data[:,dimensions[1]].max()+0.5))
    plt.show()



"""
#Takes as input a nxm numpy array,the numbers of cluster k and a fitted kmeans cluster object. KMeans object and creates a plot showing the data points with respect to the first two dimensions and colored depending on the cluster they are belonging to. Cluster centers are shown as black crosses. 
def kmeans_cluster_vis(input_data,k,kmeans_obj):
    colour = ['b','g','r','c','m','y','w']
    marker = ['o','D','v','^','s','p','*','H','+']
    
    #Checking if the selected k is within the predefined margins
    if k > (len(colour)*len(marker)):
        print("You chose a too large k. Maximum k is %s"%(len(colour)*len(marker)))
        return None
    else:
        cluster_colour = []
        for i in marker:
            for j in colour:
                cluster_colour.append(i+j)
    
    
    #Plotting data points
    for i in range(len(input_data)):
        x = input_data[i][0]
        y = input_data[i][1]
        plt.plot(x, y, cluster_colour[kmeans_obj.labels_[i]])
    
    #Plotting cluster centers
    for i in range(len(kmeans_obj.cluster_centers_)):
        x = kmeans_obj.cluster_centers_[i][0]
        y = kmeans_obj.cluster_centers_[i][1]
        plt.plot(x, y, 'kx')
        plt.text(x * (1 + 0.01), y * (1 + 0.01) , "Cluster %s"%str(i), fontsize=12)
        
    
    plt.xlim((input_data.min()-0.5, input_data.max()+0.5))
    plt.ylim((input_data.min()-0.5, input_data.max()+0.5))
    plt.show()
"""