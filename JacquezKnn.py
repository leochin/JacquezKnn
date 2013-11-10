#-------------------------------------------------------------------------------
# Name:        Space time interation test - Jacquez Knn
# Purpose:
#
# Author:      Liang-Huan Chin
#
# Created:     24/08/2013
# Copyright:   (c) Leo Chin 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
'''
This script is trying to implement Jacquez k Nearest Neighbor test
G. Jacquez. 1996. A k nearest neighbour test for space-time interaction.
            Statistics in Medicine, 15(18):1935-1949.
Two python packages are required in this script:
    1. numpy: https://pypi.python.org/pypi/numpy
    2. dbfpy: http://dbfpy.sourceforge.net/
'''

import numpy as NUM
from dbfpy import dbf as DBF

def NearestNeighbor(dict_table, kValue = 100):
    '''
    This function performs nearest neighbor search for both space and time

    input parameter: dict_table, kValue
        dict_table: an dictionary like array contains
                    ['point_number', 'X_coord', 'Y_coord', 'Time']
        kValue: int, determines how many nearest neighbors are searched

    output parameter: none, generates two csv files corresponding to space and time nearest neighbors
    no return values
    '''
    allPnt_x = NUM.array(dict_table[:,1], float)
    allPnt_y = NUM.array(dict_table[:,2], float)
    allPnt_t = NUM.array(dict_table[:,3], float)

    outputSNN = open("sNN.csv", 'w')
    outputTNN = open("tNN.csv", 'w')
    for i in xrange(arraySize):
        targetPnt_x = NUM.array(NUM.repeat(float(dict_table[i,1]), arraySize))
        targetPnt_y = NUM.array(NUM.repeat(float(dict_table[i,2]), arraySize))
        targetPnt_t = NUM.array(NUM.repeat(float(dict_table[i,3]), arraySize))

        distance = NUM.power((targetPnt_x - allPnt_x),2) + NUM.power((targetPnt_y - allPnt_y),2)
        time = NUM.abs(targetPnt_t - allPnt_t)
        outputSNN.write("Point "+str(i)+"'s nn:,")
        sNN = NUM.argsort(distance)
        for k in xrange(kValue):
            outputSNN.write(str(sNN[k])+",")
        outputSNN.write("\n")
        outputTNN.write("Point "+str(i)+"'s nn:,")
        tNN = NUM.argsort(time)
        for k in xrange(kValue):
            outputTNN.write(str(tNN[k])+",")
        outputTNN.write("\n")

    outputSNN.close()
    outputTNN.close()

def shuffle(kValue = 100):
    '''
    This function generates expected mean and variance for null hypothesis
        by shuffling X_coord, Y_coord, and Time ramdomly
        and carrying out S-T nearest neighbor search
    input parameter: kValue
        kValue: int, determines how many nearest neighbors are searched
    output parameter: none, call NearestNeighbor() and set shuffled X, Y, T as input parameter
                      generates two csv files for space and time NN
    '''
    original = dict_table[:,1]
    NUM.random.shuffle(original)
    dict_table[:,1] = original

    original = dict_table[:,2]
    NUM.random.shuffle(original)
    dict_table[:,2] = original

    original = dict_table[:,3]
    NUM.random.shuffle(original)
    dict_table[:,3] = original
    NearestNeighbor(dict_table, kValue = 100)

def cumulativeJ(sNN, tNN, kValue = 10):
    '''
    This function carries out Jacquez k NN test (cumulative sum of J)
    input parameter: sNN, tNN, kValue
        sNN: space nearest neighbor csv file path
        tNN: time nearest neighbor csv file path
        kValue: int, determines how many nearest neighbors are searched
    output parameter: an array containing a cumulative sum of Jacquez's J
    '''
    temp = 0
    cumJ = []
    for k in xrange(1,kValue+1):
        temp = temp + jacquezKnn(sNN, tNN, k)
        cumJ.append(temp)
    return cumJ

def jacquezKnn(sNN, tNN, k = 1):
    '''
    This function carries out Jacquez k NN test (specific J)
    input parameter: sNN, tNN, kValue
        sNN: space nearest neighbor csv file path
        tNN: time nearest neighbor csv file path
        kValue: int, determines how many nearest neighbors are searched
    output parameter: int, Jacquez's J at specific value of k
    '''
    sNNinput = open(sNN, 'r')
    tNNinput = open(tNN, 'r')
    J = 0
    for row in sNNinput:
        tokenS = NUM.array(row.split(','))
        tokenT = NUM.array(tNNinput.readline().split(','))
        J =  J + NUM.sum(NUM.in1d(tokenS[1:k+1], tokenT[1:k+1]))
    return J

if __name__ == '__main__':
    # this script read *.dbf as input file
    # the dbf file requires X_coord field, Y_coord field, and time_field
    path = "Infection_inRes_noWorking.dbf"
    inputFile = DBF.Dbf(path, readOnly=True)
    arraySize = len(inputFile)
    # create a dictionary-like numpy array [fid, X, Y, T]
    dict_table = NUM.array([('{0}'.format(inputFile[i][0]),
                  '{0}'.format(inputFile[i][3]),
                  '{0}'.format(inputFile[i][4]),
                  '{0}'.format(inputFile[i][5])) for i in xrange(arraySize)])
    inputFile.close()
    # call NearestNeighbor() to perform spatial and temporal nearest neighbor
    # due to limitation of memory, need to output to two csv file (space NN and time NN)
    NearestNeighbor(dict_table,50)
    # call shuffle() generate expected mean and variance
    # be careful for the duplicated name of csv file
    shuffle(10)
    # set S-T nearest neighbor csv file path
    sNN = "sNN.csv"
    tNN = "tNN.csv"
    # Jacquez's J value for k nearest neighbor
    Jvalue =  jacquezKnn(sNN, tNN, 10)
    # Cumulative sum of Jacquez's J value
    cumulative_J = cumulativeJ(sNN, tNN, kValue = 10)

##    print Jvalue
##    print cumulative_J





