#!/usr/bin/python -d

# funcENVIGEOMODELproc.py
# kingdon 2013-02-14

# ####################################################
# import modules
import os, csv, sys, argparse, re, string, itertools
import numpy as np
from numpy import *
np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)
from funcOSGEOproc import *

# ####################################################
# PARSE ENVI GEOMETRIC MODEL FILE
# accepts: inFN
# returns: xModelCoefList, yModelCoefList, polyDeg
def parseEnviModelFile(inFN):
    inFO = open(inFN,'r')
    polyCoefList = []
    for line in inFO:  # read all the lines not starting with ';' into list
        if not line.startswith(';'):  # split non-comment lines to get polynomial coefficients
            lineList = line.split()
            polyCoefList.append(lineList)
        if 'egree =' in line:  # get polynomial degree
            polyDeg = re.search(r'.*=(.*)',line)
            polyDeg = polyDeg.group(1)
            polyDeg = int(polyDeg)
            if 'ULX' in line:  # model might contain information about the base image file (see wktproj2coeffile.py)
                biULX = float(re.sub("[^\d.]+","",line))  # remove all non-numeric chars, keeping '.'
                #print biULX
            if 'ULY' in line:
                biULY = float(re.sub("[^\d.]+","",line))  # remove all non-numeric chars, keeping '.'
                #print biULY
            if 'PIXX' in line:
                biPixX = float(re.sub("[^\d.]+","",line))  # remove all non-numeric chars, keeping '.'
                #print biPixX
            if 'PIXY' in line:
                biPixY = re.sub("[^\d.]+","",line)  # remove all non-numeric chars, keeping '.'
                if '-' in line: biPixY = float('-'+biPixY)
                #print biPixY
            if 'PROJCS' in line:
                biWktPrj = line.strip(';').rstrip() # strip ; at start of line, then remove line returns
                #print biWktPrj

    #print polyCoefListl; print polyDeg
    # ####################################################
    # get length of half of the list; this will allow for separating x and y coefficients
    # flatten the list, flip strings to float, reshape array to match original input
    coefListLenHalf=len(polyCoefList)/2  
    xCoefList = ravel(polyCoefList[0:coefListLenHalf])  
    xCoefList = array([float(x) for x in xCoefList])  
    xCoefList = xCoefList.reshape(polyDeg+1,polyDeg+1)  
    yCoefList = ravel(polyCoefList[coefListLenHalf:coefListLenHalf*2])
    yCoefList = array([float(y) for y in yCoefList])
    yCoefList = yCoefList.reshape(polyDeg+1,polyDeg+1)
    # ####################################################
    xModelCoefList = []
    yModelCoefList = []
    # ####################################################
    # walk through N x N array to get 'real' polynomial coefficients
    # 'real' means exponets can't add to more than polynomial degree
    # 1st degree has 4 coefficients, 2nd has 6, 3rd has 10, ...
    for i in range(len(xCoefList)): 
        for j in range(len(xCoefList)):
            #xModelCoefList.append(xCoefList[j][i])
            if polyDeg == 1 and ((j+i) <= polyDeg or (j+i) == polyDeg + 1):
                #print i,j,xCoefList[i][j]
                xModelCoefList.append(xCoefList[i][j])
                #xModelCoefList.append(xCoefList[i][j])
            elif (j+i) <= polyDeg:
                #print j,i,xCoefList[j][i]
                xModelCoefList.append(xCoefList[j][i])
                #xModelCoefList.append(xCoefList[i][j])
    for i in range(len(yCoefList)): 
        for j in range(len(yCoefList)):
            #yModelCoefList.append(yCoefList[j][i])
            if polyDeg == 1 and ((j+i) <= polyDeg or (j+i) == polyDeg + 1):
                #print j,i,yCoefList[j][i]
                yModelCoefList.append(yCoefList[i][j])
                #yModelCoefList.append(yCoefList[i][j])
            elif (j+i) <= polyDeg:
                #print j,i,yCoefList[j][i]
                yModelCoefList.append(yCoefList[j][i])
                #yModelCoefList.append(yCoefList[i][j])
    #xModelCoefList = [-5419.20767835,1.77028151155,0.312001649093,0]
    #yModelCoefList = [13910.8252957,0.307746495151,-1.77066840229,0]
    #print xModelCoefList; print yModelCoefList
    # ####################################################
    #projectionFileTxt = open(projectionFile,'r').read().rstrip()  # read contents of .prj into a single string
    inFO.close()
    if 'PROJCS' in open(inFN,'r').read():
        biWktSRS=getWktInfo(biWktPrj)
        return (xModelCoefList, yModelCoefList, polyDeg, biULX, biULY, biPixX, biPixY, biWktSRS)
    else:
        return (xModelCoefList, yModelCoefList, polyDeg)
    #return (xModelCoefList, yModelCoefList, polyDeg, biULX, biULY, biPixX, biPixY, biWktSRS)

# ####################################################
# TRANSFORM XY POSITION
# accepts: ptXY, xCoefficientList, yCoefficientList
# returns: ptXYtrans
def polynomialXYtrans(ptXY,xCoefficientList,yCoefficientList):
    xCL=xCoefficientList; yCL=yCoefficientList
    numCoef = len(xCL)
    # hopefully all the if statements can be replaced by using a Hofman alogorithm loop
    # http://en.wikipedia.org/wiki/Horner's_method
    # http://rosettacode.org/wiki/Horner's_rule_for_polynomial_evaluation
    if len(xCL) == 4:  # RST or 1st degree polynomial
        ptXYtrans=[]
        #a1=xCL[0];a2=xCL[1];a3=xCL[2];a4=xCL[3]
        #b1=yCL[0];b2=yCL[1];b3=yCL[2];b4=yCL[3]
        # DROPPED -1 and +1 ADJUSTMENT FOR 0-BASED IMAGE ARRAY
        # this is because ENVI's GCPs are collected from 1-based image arrays
        # but IDL's prediction model for the new point location is based on  0-based image array
        ptXYtrans.append(xCL[0] + (xCL[1] * (ptXY[1])) + (xCL[2] * (ptXY[0])) + (xCL[3] * (ptXY[0]) * (ptXY[1])) )
        ptXYtrans.append(yCL[0] + (yCL[1] * (ptXY[1])) + (yCL[2] * (ptXY[0])) + (yCL[3] * (ptXY[0]) * (ptXY[1])) )

        # USE THIS MODEL FORM IF PREDICTING GCP POSITIONS THAT WERE USED IN THE CREATION OF THE ORIGINAL GEOTRANSFORM MODEL IN ENVI
        #ptXYtrans.append(xCL[0] + (xCL[1] * (ptXY[1]-1)) + (xCL[2] * (ptXY[0]-1)) + (xCL[3] * (ptXY[0]-1) * (ptXY[1]-1)) + 1)
        #ptXYtrans.append(yCL[0] + (yCL[1] * (ptXY[1]-1)) + (yCL[2] * (ptXY[0]-1)) + (yCL[3] * (ptXY[0]-1) * (ptXY[1]-1)) + 1) 

    if len(xCL) == 6:
        print("\n2nd ORDER POLYNOMIAL NOT YET READY!!\n\n")
        sys.exit(0)
        '''
        perm=itertools.permutations(range(9))
        #a1=xCL[0];a4=xCL[1];a6=xCL[2];a2=xCL[3];a5=xCL[4];a3=xCL[5]
        a1=xCL[0];a2=xCL[1];a3=xCL[2];a4=xCL[3];a5=xCL[4];a6=xCL[5]
        b1=yCL[0];b2=yCL[1];b3=yCL[2];b4=yCL[3];b5=yCL[4];b6=yCL[5]
        #[toX.append( a1 + (a2 * (pt[0]-1)) + (a3 * ((pt[0]-1)**2)) + (a4 * (pt[1]-1)) + (a5 * (pt[1]-1) * (pt[0]-1)) + (a6 * ((pt[1]-1)**2)) + 1 ) for pt in frPT]
        #[toY.append( b1 + (b2 * (pt[0]-1)) + (b3 * (pt[0]-1)**2) + (b4 * (pt[1]-1)) + (b5 * (pt[1]-1) * (pt[0]-1)) + (b6 * (pt[1]-1)**2)  + 1 ) for pt in frPT]
        #[toX.append( a1 + (a2 * (pt[1]-1)) + (a3 * ((pt[1]-1)**2)) + (a4 * (pt[0]-1)) + (a5 * (pt[0]-1) * (pt[1]-1)) + (a6 * ((pt[0]-1)**2)) + 1 ) for pt in frPT]
        #[toY.append( b1 + (b2 * (pt[1]-1)) + (b3 * (pt[1]-1)**2) + (b4 * (pt[0]-1)) + (b5 * (pt[0]-1) * (pt[1]-1)) + (b6 * (pt[0]-1)**2)  + 1 ) for pt in frPT]

        for p in perm:
            #a1=xCL[p[0]];a2=xCL[p[1]];a3=xCL[p[2]];a4=xCL[p[3]];a5=xCL[p[4]];a6=xCL[p[5]]
            #b1=yCL[p[0]];b2=yCL[p[1]];b3=yCL[p[2]];b4=yCL[p[3]];b5=yCL[p[4]];b6=yCL[p[5]]
            a1=xCL[p[0]];a2=xCL[p[1]];a3=xCL[p[2]];a4=xCL[p[3]];a5=xCL[p[4]];a6=xCL[p[5]];a7=xCL[p[6]];a8=xCL[p[7]];a9=xCL[p[8]]
            b1=yCL[p[0]];b2=yCL[p[1]];b3=yCL[p[2]];b4=yCL[p[3]];b5=yCL[p[4]];b6=yCL[p[5]];b7=yCL[p[6]];b8=yCL[p[7]];b9=yCL[p[8]]
            [toX.append( a1 + (a2 * (pt[0]-1)) + (a3 * ((pt[0]-1)**2)) + (a4 * (pt[1]-1)) + (a5 * (pt[1]-1) * (pt[0]-1)) + (a6 * ((pt[1]-1)**2)) + 1 ) for pt in frPT]
            [toY.append( b1 + (b2 * (pt[0]-1)) + (b3 * (pt[0]-1)**2) + (b4 * (pt[1]-1)) + (b5 * (pt[1]-1) * (pt[0]-1)) + (b6 * (pt[1]-1)**2)  + 1 ) for pt in frPT]
            gcpMatrix = vstack((frX,frY,toX,toY))
            gcpMatrix = gcpMatrix.T
            #print gcpMatrix[0][2]
            #sys.stdout.write('.')
            if gcpMatrix[0][2] < 663 and gcpMatrix[0][2] > 640:
                print gcpMatrix
                print p,a1,a2,a3,a4,a5,a6,a7,a8,a9
            toX=[]
            toY=[]
        '''
    if len(xCL) == 10:
        print("\n3rd ORDER POLYNOMIAL NOT YET READY!!\n\n")
        sys.exit(0)
    # 1st order: x = a1 + (a2 * (y-1)) + (a3 * (x-1)) + (a4 * (x-1) * (y-1))  + 1
    # a1 a3
    # a2 a4
    # 2nd order: x = a1 + (a2 * (y-1)) + (a3 * (y-1)^2) + (a4 * (x-1)) + (a5 * (x-1) * (y-1)) + (a6 * (x-1)^2)  + 1
    # 3nd order: x = a1 + (a2 * (y-1)) + (a3 * (y-1)^2) + (a4 * (y-1)^3) + (a5 * (x-1)) + (a6 * (x-1) * (y-1)) + 
    #                  (a7 * (x-1) * (y-1)^2) + (a8 * (x-1)^2) + (a9 * (x-1)^2 * (y-1)) + (a10 * (x-1)^3)  + 1
    #print ptXY; print ptXYtrans
    return ptXYtrans

# ####################################################
# INVERT envi GEOMETRIC MODEL
# accepts: xModelCoefList, yModelCoefList,polyDeg
# returns: xInvModelCoefList, yInvModelCoefList, polyDeg
def invertENVIgeomodel(xModelCoefList,yModelCoefList,polyDeg):
    xCL = xModelCoefList; yCL = yModelCoefList
    # create matrix of random points ... these will be used to invert the polynomial equation
    frX = random.random(1000)
    frY = random.random(1000)
    #frX = [4220,4169,3886,3957,4106]
    #frY = [4203,4162,3949,3787,4232]
    #frX = [5747.50,4879.00,5422.75,6184.00,5624.2856499]
    #frY = [2091.75,2231.25,2278.25,2048.50,2222.0324027]
    #frX = ravel(frX)+1000000
    #frY = ravel(frY)+1000000
    frPT = vstack((frX,frY)).T
    toX = []
    toY = []
    print('ENVI forward:');print(xCL); print(yCL)
    
    # call polynomialXYtrans to transform points with ENVI model coefficients
    [toX.append(polynomialXYtrans(pt,xCL,yCL)[0]) for pt in frPT]
    [toY.append(polynomialXYtrans(pt,xCL,yCL)[1]) for pt in frPT]
    gcpMatrix = vstack((frX,frY,toX,toY))
    gcpMatrix = gcpMatrix.T
    print(gcpMatrix)
    
    # calc PYTHON forward transform to compare with ENVI's coefficients
    onesArr=ones(len(frX))
    frX=ravel(frX)-1  # need to -1 for all terms for creation of models below
    frY=ravel(frY)-1
    toX=ravel(toX)-1
    toY=ravel(toY)-1
    #print toX,toY
    toXxtoY=(toX-1)*(toY-1)
    frXxfrY=(frX-1)*(frY-1)
    if len(xCL) == 4 and xCL[-1] == 0:  # if coefficient list's last element is 0 then ignore frXxfrY in model
        modelCoefMtx = vstack((onesArr,frX,frY)).T
        #print modelCoefMtx
        xSolution=np.linalg.lstsq(modelCoefMtx,toX)
        ySolution=np.linalg.lstsq(modelCoefMtx,toY)
        xCLpyforward=[xSolution[0][0],xSolution[0][2],xSolution[0][1],0]  # switch [0][1] and [0][2] to match ENVI's order
        yCLpyforward=[ySolution[0][0],ySolution[0][2],ySolution[0][1],0]
    if len(xCL) == 4 and not xCL[-1] == 0:  # if not 0 then include frXxfrY in model
        modelCoefMtx = vstack((onesArr,frX,frY,frXxfrY)).T
        xSolution=np.linalg.lstsq(modelCoefMtx,toX)
        ySolution=np.linalg.lstsq(modelCoefMtx,toY)
        xCLpyforward=[xSolution[0][0],xSolution[0][2],xSolution[0][1],xSolution[0][3]]  # switch [0][1] and [0][2] to match ENVI's order
        yCLpyforward=[ySolution[0][0],ySolution[0][2],ySolution[0][1],ySolution[0][3]]
    print('PYTHON forward:'); print(xCLpyforward); print(yCLpyforward)
    
    # calc PYTHON backward: invert model by regressing toX and toY against frX and frY, then re-calc model coefficients
    if len(xCL) == 4 and xCL[-1] == 0:  # if coefficient list's last element is 0 then ignore frXxfrY in model
        modelType='rst'
        modelCoefMtx = vstack((onesArr,toX,toY)).T
        xSolution=np.linalg.lstsq(modelCoefMtx,frX)
        ySolution=np.linalg.lstsq(modelCoefMtx,frY)
        xCLpybackward=[xSolution[0][0],xSolution[0][2],xSolution[0][1],0]  # switch [0][1] and [0][2] to match ENVI's order
        yCLpybackward=[ySolution[0][0],ySolution[0][2],ySolution[0][1],0]
    if len(xCL) == 4 and not xCL[-1] == 0:  # if not 0 then include frXxfrY in model
        modelType='py1'
        modelCoefMtx = vstack((onesArr,toX,toY,toXxtoY)).T
        xSolution=np.linalg.lstsq(modelCoefMtx,frX)
        ySolution=np.linalg.lstsq(modelCoefMtx,frY)
        xCLpybackward=[xSolution[0][0],xSolution[0][2],xSolution[0][1],xSolution[0][3]]  # switch [0][1] and [0][2] to match ENVI's order
        yCLpybackward=[ySolution[0][0],ySolution[0][2],ySolution[0][1],ySolution[0][3]]
    print('PYTHON backward:',modelType); print(xCLpybackward); print(yCLpybackward)
    
    #'''
    # apply backward coefficients to gcps to check if inversion worked
    # call polynomialXYtrans to transform point
    toX=ravel(toX)+1  # +1 to values to make them 'real' again, then do backward transform
    toY=ravel(toY)+1
    toPT = vstack((toX,toY)).T
    frX=[];frY=[]
    [frX.append(polynomialXYtrans(pt,xCLpybackward,yCLpybackward)[0]) for pt in toPT]
    [frY.append(polynomialXYtrans(pt,xCLpybackward,yCLpybackward)[1]) for pt in toPT]
    gcpMatrix = vstack((toX,toY,frX,frY))
    gcpMatrix = gcpMatrix.T
    print(gcpMatrix)
    #print gcpMatrix[range(len(gcpMatrix))]
    #'''
    
    return (xCLpybackward,yCLpybackward)
    
    
