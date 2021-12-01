#!/usr/bin/python3 -d

# imgpoint.py
# kingdon 2021-09-15

# ####################################################
# import modules
import os, os.path, sys, time, argparse, csv, re, math
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
from funcOSGEOproc import *
from funcENVIGEOMODELproc import *

# get current dir
cwd = os.getcwd()
#print 'cwd: %s' %(cwd)

# ####################################################
# MANAGE ARGUMENTS WITH argparse
parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
        #formatter_class = argparse.RawTextHelpFormatter,
    description='extract pixel values based on point locatons',
    epilog='\ndepends on the following python modules: \nosgeo (gdal, ogr, osr, gdalconst), funcOSGEOproc.py (home-rolled) \
            \nos, os.path, sys, argparse, csv, re\n\n\n')
# using ALL optional arguments just to keep things neat, though all args will be required
parser.add_argument('-p','--point-input',type=str,nargs='?',help='input point file; either .shp or .csv')
parser.add_argument('-r','--raster-input',type=str,nargs='?',help='input raster file; anything supported by GDAL')
parser.add_argument('-c','--csv-output',type=str,nargs='?',help='csv file created to hold extracted pixel values')
# 'real' optional arguments
parser.add_argument('-x','--xcoord-field',type=str,nargs='?',help='optional: point file\'s X coordinate field; only used if POINT_INPUT is a .csv')
parser.add_argument('-y','--ycoord-field',type=str,nargs='?',help='optional: point file\'s Y coordinate field; only used if POINT_INPUT is a .csv')
parser.add_argument('-f','--fid-field',type=str,nargs='?',help='optional: point file field containing UNIQUE identifier; if not specified, a field named ZID will be created in CSV_OUTPUT')
parser.add_argument('-e','--epsg-code',type=str,nargs='?',help='optional: input EPSG code for point file; .prj associated with .shp takes precedence')
#parser.add_argument('-k','--keep-all',type=str,nargs='?',default='yes',help='optional: <yes|no>; default=yes; keep all point even if they don\'t intersect raster')
parser.add_argument('-k','--keep-all',action='store_true',help='NOT WORKING YET!! optional: default=yes; keep all point even if they don\'t intersect raster')
parser.add_argument('-n','--nodata-value',type=int,nargs='?',help='NOT WORKING YET!! optional: value assigned to points that don\'t intersect RASTER_FILE')
parser.add_argument('-i','--ignore-rotation',action='store_true',help='optional: ignore rotation in ENVI header')
parser.add_argument('-t','--transformation-file',type=str,nargs='?',help='optional: output from wktproj2coeffile.py; modified ENVI polynomial coefficient transformation file')
parser.add_argument('-4','--four-nearest',action='store_true',help='optional:collect values for 4 nearest pixels to point')
args=parser.parse_args()
inPts,inRas,outCsv,xFld,yFld,fFld,epsgCode,keepAll,noData,ignoreRot,coefFile,fourNear = \
[args.point_input,args.raster_input,args.csv_output,args.xcoord_field,args.ycoord_field,args.fid_field,args.epsg_code, \
args.keep_all,args.nodata_value,args.ignore_rotation,args.transformation_file,args.four_nearest]
#print
#print args
#print
#print inPts,inRas,outCsv,xFld,yFld,fFld,epsgCode,ignoreRot,keepAll,noData,coefFile
#print ignoreRot,coefFile
print()

# error checking arguments
if None in [inPts,inRas,outCsv]:
    print('\n !!! missing one or more mandatory arguments !!! \n')
    parser.print_help()
    sys.exit(1)

# time starts
startTime=time.time()

# manipulate args that are file names
if not inPts == None:   inPts=os.path.abspath(inPts)
if not inRas == None:   inRas=os.path.abspath(inRas)
if not coefFile == None:    coefFile=os.path.abspath(coefFile)
#if not outCsv == None and not '/' in outCsv:   outCsv=cwd+'/'+outCsv
if not outCsv == None and not '/' in outCsv:
    outCsv=cwd+'/'+outCsv
else:
    outCsv=os.path.abspath(outCsv)

# is inPts a shape file or csv?
doShp=False; doCsv=False
inPtsExt = inPts[-4:].lower()
if inPtsExt == '.shp':   doShp=True; inShpPrj=inPts[:-4]+'.prj' #print 'shp: '+inPtsExt
if inPtsExt == '.csv':   doCsv=True; #print 'csv: '+inPtsExt
#print inShpPrj

# is an inverse polynomial transform required?
doInvPolyTransform=False
if os.path.isfile(str(coefFile)): doInvPolyTransform=True # NEED TO WRAP FILENAME IN STR() !!!

# error checking arguments
if not os.path.isfile(inPts) or not os.path.isfile(inRas):
    print('\n !!! either '+inPts+' or '+inRas+' does not exist !!!\n')
    parser.print_help()
    sys.exit(1)
if not os.path.exists(os.path.dirname(outCsv)):
    print('\n !!! path to '+outCsv+' does not exist !!!\n')
    parser.print_help()
    sys.exit(1)
if not os.path.exists(os.path.dirname(outCsv)):
    print('\n !!! path to '+outCsv+' does not exist !!!\n')
    parser.print_help()
    sys.exit(1)
if doCsv and None in [xFld,yFld,epsgCode]:
    print('\n !!! POINT_INPUT is a CSV file !!!\n !!! you must specify XCOORD_FIELD YCOORD_FIELD EPSG_CODE !!!\n')
    parser.print_help()
    sys.exit(1)
if doShp and not os.path.exists(inShpPrj) and not epsgCode:
    print('\n !!! POINT_INPUT is a SHP file !!!\n !!! but there is no associated .prj file !!! \n you must specify an EPSG_CODE or create a .prj !!!\n')
    parser.print_help()
    sys.exit(1)
if doShp and os.path.exists(inShpPrj) and epsgCode:
    print('\n WARNING: using .prj file instead of specified EPSG_CODE')


# ####################################################
# GET INFORMATION ABOUT inRas
#getRasterInfo(inRas) # call to getRasterInfo
[inRasCols, inRasRows, inRasMinX, inRasMaxY, inRasMaxX, inRasMinY, \
inRasCtrX, inRasCtrY, inRasPixX, inRasPixY, inRasNumBands, inRasRotAng, inRasGTRotAng,\
inRasGeoTransform, inRasSRS, inRasPrjWkt, inRasFileType,inRasDataType] = getRasterInfo(inRas)
inRasBandFieldNames = [''.join(["band_",str(i+1)])for i in range(inRasNumBands)]
inRasBandFieldNames = ','.join(inRasBandFieldNames)
inRasRotAngRad = -1*(float(inRasRotAng)*math.pi/180)
#inRasRotAngRad = 0
#print inRasNumBands, inRasBandFieldNames
#print inRasMinX, inRasMaxY, inRasMaxX, inRasMinY, inRasCtrX, inRasCtrY, inRasRotAng
#print 'inRasSRS: ' + str(inRasSRS)
#print isinstance(inRasMinX,int); print isinstance(inRasMinX,float)

# inRas has to be read again because layer can't be passed back from getRasterInfo function
inRasDriver = gdal.AllRegister
# open the file
inRasDS = gdal.Open(inRas,GA_ReadOnly)

# ####################################################
# GET INFORMATION ABOUT inPts, DIFFERENT IF shp OR csv 
if doShp:
    #inShpPrj = (inPts)[:-4]+'.prj'
    ([inPtsExtent,inPtsNumPts,inPtsFields,inPtsSRS]) = getShapeFileInfo(inPts,inShpPrj)
    if not  os.path.exists(inShpPrj):   ([inPtsSRS])=getEPSGInfo(int(epsgCode)) # get SRS for EPSG code
    # shape file has to be read again because layer isn't being passed back from getShapeFileInfo function
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    inPtsDS = shpDriver.Open(inPts)
    inPtsLayer = inPtsDS.GetLayer()
    inPtsLayerDefn = inPtsLayer.GetLayerDefn()
if doCsv:
    ([inPtsSRS])=getEPSGInfo(int(epsgCode)) # get SRS for EPSG code
    inPtsCSVFile = open(inPts,'rb')
    inPtsDS = csv.reader(inPtsCSVFile)
    inPtsFields = next(inPtsDS) # read first line from incoming CSV
    inPtsLayer = inPtsDS
    xFldIdx = inPtsFields.index(xFld) # get index of X field
    yFldIdx = inPtsFields.index(yFld) # get index of Y field
    #print 'x/y idx: ' + str(xFldIdx) + ',' + str(yFldIdx)
    # get number of points (lines) in file
    i = -1
    with open(inPts) as f:
            for i, l in enumerate(f):
                    pass
    inPtsNumPts = i


#print 'inPtsSRS: ' + str(inPtsSRS)+'\n'
#print 'inPtsFields: ' + str(inPtsFields)

############################################################
# GET fFld FROM inPtsFields
try:
    fFldIdx = inPtsFields.index(fFld)
except ValueError:
    fFldIdx = -1
imgRowA='9999'
if fFld == None:
    outCsvFields = 'ZID,imgX,imgY,imgXrot,imgYrot,imgCol,imgRow,imgRowA,ptX,ptY,'+str(inRasBandFieldNames)
elif fFldIdx > -1:
    outCsvFields = 'ZID,'+fFld+',imgX,imgY,imgXrot,imgYrot,imgCol,imgRow,imgRowA,ptX,ptY,'+str(inRasBandFieldNames)
elif fFldIdx == -1:
    print('\n !!! '+fFld+' is not a field in '+inPts+' !!!\n')
    parser.print_help()
    sys.exit(1)
#print outCsvFields 

# ####################################################
# GET INFORMATION ABOUT coefFile's SRS
# PUT ALL OF THIS INTO A FUNCTION IN funcENVIGEOMODELproc.py
if doInvPolyTransform == True:
    coefFileInfo=parseEnviModelFile(coefFile)
    #print coefFileInfo
    [xModelCoefList, yModelCoefList, polyDeg, biULX, biULY, biPixX, biPixY, biWktSRS] = coefFileInfo
    coefFileSRS = biWktSRS

# ####################################################
# DO COORDINATE TRANSFORMATIONS 
coordTrans = osr.CoordinateTransformation(inPtsSRS,inRasSRS) # create geographic transformation for points
if doInvPolyTransform == True: # create transform from inPtsSRS to SRS based on WKT in ENVI coefFile transformation file
    coordTrans = osr.CoordinateTransformation(inPtsSRS,coefFileSRS)
    
    
# ####################################################
# CREATE OUTPUT CSV FILE, GET READY TO STEP THROUGH POINTS 
outCsvFH = open(outCsv,'w')  # create output csv file
outCsvFH.write(outCsvFields+'\n')  # write first line to output csv file

fidVal = 0
ptsInRasExtents = 0
if doShp:   pt=inPtsLayer.GetNextFeature()
if doCsv:   pt=next(inPtsLayer)
#print pt

# set raster extents to file column/row for polynomial transform
if doInvPolyTransform == True:
    inRasMinX = 0
    inRasMinY = 0
    inRasMaxX = inRasCols
    inRasMaxY = inRasRows 
#print inRasRows
#print biULX,biULY
#print inRasMinX,inRasMaxX
#print inRasMinY,inRasMaxY
#print
############################################################
# LOOP THROUGH POINT FEATURES (shp) / FILE LINES (csv)
while pt:
    #print 'PT'
    if doShp:
        ptgeom=pt.GetGeometryRef()
        xOriginal=ptgeom.GetX()
        yOriginal=ptgeom.GetY()
    if doCsv:
        xOriginal=pt[xFldIdx]
        yOriginal=pt[yFldIdx]
        ptgeom=ogr.Geometry(ogr.wkbPoint)
        ptgeom.AddPoint(float(xOriginal),float(yOriginal))

    ptgeom.Transform(coordTrans) # 'reproject' the point
    imgX=ptgeom.GetX()
    imgY=ptgeom.GetY()
            
    # if rotAng <> 0, then rotate point around raster origin
    #if not inRasRotAng == 0 and not ignoreRot == 'yes' and doInvPolyTransform == False:
    if not inRasRotAng == 0 and not ignoreRot and doInvPolyTransform == False:
        #print 'ROTATION'
        # shift point so that origin is at raster origin 
        imgXshift=float(imgX) - float(inRasMinX)
        imgYshift=float(imgY) - float(inRasMaxY)
        # do the rotation
        imgXrot=imgXshift*math.cos(float(inRasRotAngRad)) - imgYshift*math.sin(float(inRasRotAngRad))
        imgYrot=imgXshift*math.sin(inRasRotAngRad) + imgYshift*math.cos(inRasRotAngRad)
        # shift the points back to their original origin
        imgXrot = imgXrot + float(inRasMinX)
        imgYrot = imgYrot + float(inRasMaxY)
        #print '--',imgXrot, imgYrot
    elif doInvPolyTransform == True:
        #print 'POLYNOMIAL TRANFORMATION'
        # biPixY is -ve, so it gets added
        # shift pixels to match base image from ENVI model file
        imgXrot=((imgX - biULX) / biPixX)
        imgYrot=((imgY - biULY) / biPixY)
        # do polynomial transform to move point from map space to file space
        [imgXrot,imgYrot] = polynomialXYtrans([imgXrot,imgYrot],xModelCoefList,yModelCoefList)
        #print imgXrot,imgYrot
    else:
        #print 'NO ROTATION, NO POLYNOMIAL TRANFORMATION'
        imgXrot=imgX
        imgYrot=imgY
    
    #if fFld > -1:   
    if fFld is not None:   
        if doShp:   fFldVal = pt.GetFieldAsString(fFld)  # need to get field value before advancing the pt
        if doCsv:   fFldVal = str(pt[fFldIdx])  # need to get field value before advancing the pt
    fidVal = fidVal + 1
    
        
    if fourNear == False:  # just being explicit here
        #print imgX, imgY, imgXrot, imgYrot
        print(imgXrot, imgYrot, inRasMinX, inRasPixX)
        if (not math.isinf(imgXrot)) and (not math.isinf(imgYrot)):
            if doInvPolyTransform == False:
                imgCol=int((imgXrot - inRasMinX) / inRasPixX)
                imgRow=int((imgYrot - inRasMaxY) / inRasPixY)
                #print 'no poly:',imgCol,imgRow
            elif doInvPolyTransform == True:  
                #imgCol = int(imgXrot-0.5)
                #imgRow = int(imgYrot+0.5)
                #imgRowA = inRasRows-imgRow
                # when compairing corrected and non-correct points over 6 images (~150 points), -0.5 pixel shift minimized mean bias
                imgCol = int(imgXrot-0.5)
                imgRow = int(imgYrot-0.5)
                imgXrot = imgXrot - 0.5  # this adjusts the value to better match the pixel that is actually sampled; must happen after imgCol adjustment
                imgRowA = inRasRows-imgYrot + 1.0 - 0.5  # 1.0 accounts for starting row; 0.5 accounts for pixel shift (as for imgCol, imgRow)
                # imgXrot, imgYrot, and imgRowA are only used for displaying points later; they're not used for extracting raster pixel value
        # if FID field has been specified, get the value and write it to the output CSV
        ################# ADDED imgRowA HERE !!!! CAN REMOVE IT LATER, AND ABOVE WHERE CSV COLUMNS ARE NAMED!!!########################
        # only extract pixel values if point falls within the raster 
        #print fFldVal,imgXrot,imgYrot; print inRasMinX, inRasMaxX; print inRasMinY, inRasMaxY
        if inRasMinX < imgXrot < inRasMaxX and inRasMinY < imgYrot < inRasMaxY:
            #if fFld > -1:
            if fFld is not None:   
                s = ','.join([str(x) for x in [fidVal,fFldVal,imgX,imgY,imgXrot,imgYrot,imgCol,imgRow,imgRowA,xOriginal,yOriginal]])
            else:
                s = ','.join([str(x) for x in [fidVal,imgX,imgY,imgXrot,imgYrot,imgCol,imgRow,imgRowA,xOriginal,yOriginal]])
    
            # do int on imgCol and imgRow to slam point to pixel's upper left corner
            imgCol = int((imgCol))
            imgRow = int((imgRow))
            #print 'readasarray',imgCol, imgRow
            for b in range(inRasNumBands):
                band = inRasDS.GetRasterBand(b+1)  # 1-based array index
                pixVal = band.ReadAsArray(imgCol,imgRow,1,1) # returns 2-D array
                pixVal = pixVal[0,0] # get value 0,0 from that 2-D array
                s = s + ',' + str(pixVal)
                #print fFldVal,ptsInRasExtents
            ptsInRasExtents = ptsInRasExtents + 1
            #print s
            outCsvFH.write(s+'\n')
    
    # FOUR NEAREST PIXELS
    elif fourNear == True:
        #print fFldVal,imgX, imgY
        # find coords at NE SE SW NW or 45 135 225 315
        for dirAng in (45,135,225,315):
            if doInvPolyTransform == False:
                imgXrotOriginal=imgXrot
                imgYrotOriginal=imgYrot
                imgXrot = imgXrotOriginal
                imgYrot = imgYrotOriginal
                pixSize=inRasPixX
                #pixSize = pixSize * 0.5
                #print pixSize
                aRad = -1*(float(dirAng-90)*math.pi/180)  # subtracting 90 from dirAng makes 0 == North ?
                pixHypHalf = (((pixSize**2)+(pixSize**2))**0.5) / 2 # use pixel size to calc 1/2 hypotenuse
                #print pixHypHalf
                ########## ADD URL HERE ###############
                # 0 degrees = NORTH
                #imgXrot=imgXrot+pixHypHalf*math.cos(aRad)-pixHypHalf*math.sin(aRad)
                #imgYrot=imgYrot+pixHypHalf*math.sin(aRad)+pixHypHalf*math.cos(aRad)
                imgXrot=imgXrot+pixHypHalf*math.cos(aRad)
                imgYrot=imgYrot+pixHypHalf*math.sin(aRad)
                #print 'post rot:',imgXrot,imgYrot
                imgCol=int((imgXrot - inRasMinX) / inRasPixX)
                imgRow=int((imgYrot - inRasMaxY) / inRasPixY)
                #print 'no poly:',imgCol,imgRow
            elif doInvPolyTransform == True:  
                # shift pixels to match base image from ENVI model file
                imgXrotOriginal=((imgX - biULX) / biPixX)
                imgYrotOriginal=((imgY - biULY) / biPixY)
                # store point location; will reset this to do each rotation in loop
                imgXrot = imgXrotOriginal
                imgYrot = imgYrotOriginal
                #print 'pre  rot:',imgXrot,imgYrot
                # PIXEL SIZE USED HERE IS NOT biPixX; HYPOTENUSE NEEDS TO BE BASED ON inRasPixX, BUT TRANSFORMED
                #pixSize=1*0.5
                pixSize = float(inRasPixX / biPixX)
                #pixSize = pixSize*0.5
                aRad = (float(dirAng-90)*math.pi/180)  # adding 180 to dirAng makes 0 == North ?
                pixHypHalf = (((pixSize**2)+(pixSize**2))**0.5) / 2 # use pixel size to calc 1/2 hypotenuse
                #print
                #print 'inRasPixX',inRasPixX
                #print 'biPixX', biPixX
                #print 'pixSize',pixSize
                #print 'pixHypHalf',pixHypHalf
                ########## ADD URL HERE ###############
                # 0 degrees = NORTH
                imgXrot=imgXrot+pixHypHalf*math.cos(aRad)
                imgYrot=imgYrot+pixHypHalf*math.sin(aRad)
                #imgXrot=imgXrot+pixSize*math.cos(aRad)
                #imgYrot=imgYrot+pixSize*math.sin(aRad)
                #imgXrot=imgXrot+pixHypHalf*math.cos(aRad)+pixHypHalf*math.sin(aRad)
                #imgYrot=imgYrot+pixHypHalf*math.sin(aRad)-pixHypHalf*math.cos(aRad)
                #imgXrot=imgXrot+pixSize*math.cos(aRad)+pixSize*math.sin(aRad)
                #imgYrot=imgYrot+pixSize*math.sin(aRad)-pixSize*math.cos(aRad)
                [imgXrot,imgYrot] = polynomialXYtrans([imgXrot,imgYrot],xModelCoefList,yModelCoefList)
                #print 'post rot:',imgXrot,imgYrot
                # when compairing corrected and non-correct points over 6 images (~150 points), -0.5 pixel shift minimized mean bias
                imgCol = int(imgXrot-0.5)
                imgRow = int(imgYrot-0.5)
                imgXrot = imgXrot - 0.5  # this adjusts the value to better match the pixel that is actually sampled; must happen after imgCol adjustment
                imgRowA = inRasRows-imgYrot + 1.0 - 0.5  # 1.0 accounts for starting row; 0.5 accounts for pixel shift (as for imgCol, imgRow)
                # imgXrot, imgYrot, and imgRowA are only used for displaying points later; they're not used for extracting raster pixel value
            # if FID field has been specified, get the value and write it to the output CSV
            ################# ADDED imgRowA HERE !!!! CAN REMOVE IT LATER, AND ABOVE WHERE CSV COLUMNS ARE NAMED!!!########################
            # only extract pixel values if point falls within the raster 
            #print fFldVal,imgXrot,imgYrot; print inRasMinX, inRasMaxX; print inRasMinY, inRasMaxY
            if inRasMinX < imgXrot < inRasMaxX and inRasMinY < imgYrot < inRasMaxY:
                if fFld > -1:
                    fFldValPfx = fFldVal
                    fFldVal=fFldVal+"_"+str(dirAng)
                    s = ','.join([str(x) for x in [fidVal,fFldVal,imgX,imgY,imgXrot,imgYrot,imgCol,imgRow,imgRowA,xOriginal,yOriginal]])
                    fFldVal = fFldValPfx
                else:
                    fidValPfx = fidVal
                    s = ','.join([str(x) for x in [fidVal,imgX,imgY,imgXrot,imgYrot,imgCol,imgRow,imgRowA,xOriginal,yOriginal]])
                    fidVal = fidValPfx
                # do int on imgCol and imgRow to slam point to pixel's upper left corner
                imgCol = int((imgCol))
                imgRow = int((imgRow))
                #print 'readasarray',imgCol, imgRow
                for b in range(inRasNumBands):
                    band = inRasDS.GetRasterBand(b+1)  # 1-based array index
                    pixVal = band.ReadAsArray(imgCol,imgRow,1,1) # returns 2-D array
                    pixVal = pixVal[0,0] # get value 0,0 from that 2-D array
                    s = s + ',' + str(pixVal)
                ptsInRasExtents = ptsInRasExtents + 0.25

                #print s
                outCsvFH.write(s+'\n')
            #print

    

    # move to next point
    #ptsInRasExtents = ptsInRasExtents + 1



    if doShp:   pt.Destroy(); pt=inPtsLayer.GetNextFeature()
    if doCsv:   
        try:
            pt=next(inPtsLayer)  # insert try statement to catch StopIteration error
        except StopIteration:
            break

if doShp:   inPtsDS.Destroy()

ptsInRasExtents = int(ptsInRasExtents)
############################################################
# END TIME
endTime = time.time()
totalTime = str(endTime - startTime) + ' seconds'
print('\n\n '+str(ptsInRasExtents)+'/'+str(inPtsNumPts)+' points from '+os.path.basename(inPts))
print('\n are within the extents of raster\n\n '+os.path.basename(inRas))
print('\n '+ 'total time: '+totalTime+'\n\n')
if keepAll == True or not noData == None:
    print('you used the -k or -n options, or both\n')
    print('please note that neither of these options are operational\n')
    print('imgpoints.py completed successfully, but ignored these options\n\n')



############################################################
# CLEAN UP
outCsvFH.close()
if doCsv:   inPtsCSVFile.close()

############################################################
#sys.exit(1)
print()
print()
############################################################













