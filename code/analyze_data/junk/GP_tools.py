## Tools for GP

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import pickle


# ML imports
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics.pairwise import euclidean_distances
from scipy.spatial import distance
from scipy import optimize, linalg
import scipy


def linear_kernel(X, X_, hypers):
    """ The linear (dot product) kernel for two inputs. """
    vp = hypers[0]
    return vp * X @ X_.T
	

def linear_kernel(X, X_, hypers):
    """ Calculate the Matern kernel between X and X_.
    Parameters:
        X (np.ndarray):
        X_ (np.ndarray)
        hypers (iterable): default is ell=1.0.
    Returns:
        K (np.ndarray)
    """
    #D = distance.cdist(X, X_)
    D = euclidean_distances(X, X_)
    D_L = D / hypers[0]

    first = (1.0 + np.sqrt(5.0) * D_L) + 5.0 * D_L ** 2 / 3.0
    second = np.exp(-np.sqrt(5.0) * D_L)

    K = first * second
    return K


def predict_GP(X_train, y_train, X_test, prams):
    """ Gaussian process regression predictions.
    Parameters:
        X_train (np.ndarray): n x d training inputs
        y_train (np.ndarray): n training observations
        X_test (np.ndarray): m x d points to predict
    Returns:
        mu (np.ndarray): m predicted means
        var (np.ndarray): m predictive variances
    """

    # Evaluate kernel on training data
    K = linear_kernel(X_train, X_train, prams[1:])

    # To invert K_y we use the Cholesky decomposition (L)
    L = np.linalg.cholesky(K + np.eye(np.shape(X_train)[0])*prams[0]**2)

    # solve for z=L^-1y
    z = linalg.solve_triangular(L, y_train, lower=True)
    alpha = linalg.solve_triangular(L.T, z, lower=False)
    K_star = linear_kernel(X_train, X_test, prams[1:])
    mu = np.matmul(K_star.T, alpha)

    # Compute the variance at the test points
    z = linalg.solve_triangular(L, K_star, lower=True)
    alpha = linalg.solve_triangular(L.T, z, lower=False)
    K_star_star = linear_kernel(X_test, X_test, prams[1:])
    v = np.diag(K_star_star) - np.dot(K_star.T, alpha)
    v = np.diag(v)
    return mu, v

def neg_log_marg_likelihood(log_prams, X, y):
    """ Calculate the negative log marginal likelihood loss.
    We pass the log hypers here because it makes the optimization
    more stable.
    Parameters:
        log_hypers (np.ndarray): natural log of the hyper-parameters
        X (np.ndarray)
        y (np.ndarray)
    Returns:
        (float) The negative log marginal likelihood.
    """

    non_log_prams = np.exp(log_prams)
    #print(non_log_prams)

    # Evaluate kernel on training data
    K = linear_kernel(X, X, non_log_prams[1:])

    # To invert K we use the Cholesky decomposition (L), because symmetric and positive  definite
    n = len(y)
    L = np.linalg.cholesky(K + np.eye(np.shape(X)[0])*non_log_prams[0])
    z = linalg.solve_triangular(L, y, lower=True)
    alpha = linalg.solve_triangular(L.T, z, lower=False) #dont know about this

    log_p_y_X = 0.5*np.matmul(y, alpha) + np.sum(np.log(np.diag(L))) + 0.5*n*np.log(2*np.pi)
    return log_p_y_X
