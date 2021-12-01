# module funcOSGEOproc

# ####################################################
# import modules
import os, os.path, sys, time, argparse, csv, re
from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
from math import *

# getRasterInfo  rF is rasterFile
def getRasterInfo(rF):
	rFDriver = gdal.AllRegister
	# open the file
	rFDS = gdal.Open(rF,GA_ReadOnly)
	# columns, rows, bands
	rFCols = rFDS.RasterXSize
	rFRows = rFDS.RasterYSize
	rFBands = rFDS.RasterCount
	rFDataType = gdal.GetDataTypeName(rFDS.GetRasterBand(1).DataType)
	# ulX, ulY, pixX, pixY
	rFGeoTransform = rFDS.GetGeoTransform()
	#print rFGeoTransform
	# get lrX, lrY
	[rFULX, rFGT1X, rFGT2X, rFULY, rFGT4Y, rFGT5Y] =  rFGeoTransform
	#rFULX = rFGeoTransform[0]; rFULY = rFGeoTransform[3]; rFPixX = rFGeoTransform[1]; rFPixY = rFGeoTransform[5]  # another way to get tranform parameters
	#print 'rotated'
	rFLRX = rFULX + rFCols*rFGT1X + rFRows*rFGT2X
	rFLRY = rFULY + rFCols*rFGT4Y + rFRows*rFGT5Y
	rFCtrX = rFULX + ((rFLRX - rFULX) / 2)
	rFCtrY = rFULY - ((rFULY - rFLRY) / 2)
	# calc rotation angle
	rFGTRotAng = degrees(atan(rFGT2X/rFGT1X))
	#print rFGTRotAng
	# calc pixel size
	if rFGTRotAng != 0:
		rFPixX = rFGT1X / cos(radians(rFGTRotAng))
		rFPixY = -1 * (rFGT4Y / sin(radians(rFGTRotAng))) 
	else:
		rFPixX = rFGT1X
		rFPixY = rFGT5Y
	# all this (above) borks if angle is 90 or 270 because cos(90) or cos(270) is 0
	# so, I do the check below and pick the bigger of the two calculated pixels sizes
	# really, there should never be a reason (i.e. no benefit) to use a 90,180,270, etc rotation on a non-north-up image
	# if 90 rotation makes it north-u then just store it this way!!
	if (float(rFGTRotAng) == 90 or float(rFGTRotAng == 270)):
		absPix = max(abs(rFPixX),abs(rFPixY))
		rFPixX = absPix
		rFPixY = -1 * absPix
	#print 'rFPixX', rFPixX; print 'rFPixY', rFPixY
	#print rFCols, rFRows, rFPixX, rFPixY
	#print rFULX, rFULY, rFLRX, rFLRY
	#rFULX = "%.3f" % rFULX
	# get raster file type
	rFDriver = rFDS.GetDriver()
	#rFFileType = rFFileType.ShortName 
	rFFileType = rFDriver.LongName 
	# check if raster is ENVI, has a .hdr, then look for rotation factor
	# this is only for rasters with rotation that IS NOT represented in geotranform parameters
	if 'ENVI' in rFFileType:
		#rFHdr =rF[:-4]+'.hdr'
		#rFHdr = os.path.basename(os.path.abspath(rF))
		#print rF
		try:	# this try/except is to catch the IOError that may arise if rF's extension is .hdr or, say, .bsq.hdr
			#rFHdr = os.path.splitext(rFHdr)[0]+'.hdr'
			#rFHdrTxt = open(rFHdr,'r').read()
			rFHdr = os.path.splitext(rF)[0]+'.hdr'
			#rFHdrTxt = open(rF,'r').read()
			rFHdrTxt = open(rFHdr,'r').read()
		except IOError:
			rFHdr = rF+'.hdr'
			rFHdrTxt = open(rFHdr,'r').read()
		#`sed -n 's/.*rotation.*=\(.*\)\}.*/\1/p' $envihdr` # this is how it would work in sed 
		#p = re.compile(r'rotation.*=(.*).*}')  # could make expression with this line
		pattern = r'rotation.*=(.*).*}' # get text between = and }, which is the rotation angle
		match = re.search(pattern,rFHdrTxt)
		if match:	# remember this for later when you want to do point rotation
			rFRotAng = match.group(1)
			#print 'rotAng: '+rotAng
		else:
			rFRotAng = 0
	else:
		rFRotAng = 0
	# spatial ref system
	rFPrjWkt = rFDS.GetProjection() # get the rF's projection
	rFSRS = osr.SpatialReference()
	rFSRS.ImportFromWkt(rFPrjWkt) # import the spatial reference system from rFPrjWkt (assumes returned proj is well-known text)
	return ([rFCols, rFRows, rFULX, rFULY, rFLRX, rFLRY, rFCtrX, rFCtrY, rFPixX, rFPixY, rFBands, rFRotAng, rFGTRotAng, rFGeoTransform,rFSRS,rFPrjWkt,rFFileType,rFDataType])

# getShapeFileInfo  sF is shapeFile
def getShapeFileInfo(sF,projectionFile):
	#print 'from doShpFunc: '+sF
	#print 'from doShpFunc: '+projectionFile
	try:
		projectionFileTxt = open(projectionFile,'r').read().rstrip()  # read contents of .prj into a single string
		# spatial ref system
		sFSRS = osr.SpatialReference()
		sFSRS.ImportFromWkt(projectionFileTxt) # import the spatial reference system from sFProj (assumes returned proj is well-known text); though it came from an ESRI .prj, I didn't user ImportFromESRI)
		#print 'sFSRS: ' + str(sFSRS)
		#projectionFileTxt.close()
	except IOError:
		sFSRS = None
		pass
	#print projectionFileTxt
	shpDriver = ogr.GetDriverByName("ESRI Shapefile")
	sFDS = shpDriver.Open(sF)
	sFLayer = sFDS.GetLayer()
	sFLayerDefn = sFLayer.GetLayerDefn()
	sFExtent = sFLayer.GetExtent()
	sFNumPts = sFLayer.GetFeatureCount()
	sFFields=[sFLayerDefn.GetFieldDefn(i).GetName() for i in range(sFLayerDefn.GetFieldCount())]
	#sFFields=','.join(sFFields)  # don't convert sFFields to string
	#print 'number of point features: ' + str(sFNumPts)
	#print 'shapefile extent:' + str(extent)
	#return ([sFExtent,sFNumPts,sFFields,sFSRS,sFLayer])
	return ([sFExtent,sFNumPts,sFFields,sFSRS])
	#sFDS.Destroy()

# getEPSGCodeSRS  
def getEPSGInfo(eCode):
	eCodeSRS = osr.SpatialReference()
	eCodeSRS.ImportFromEPSG(eCode)
	return ([eCodeSRS])

# getWktSRS  
def getWktInfo(wktString):
	wktSRS = osr.SpatialReference()
	wktSRS.ImportFromWkt(wktString)
	#print wktSRS
	#return ([wktSRS])
	return (wktSRS)


# getBlockSize
#def getBlockSize(rF):
	#x = gdal.ptrcreate('int',0,2)
	#gdal.GDALGetBlockSize(band

# calcProjWin
# accepts: ulx uly lrx lry rasMinX rasMaxY pixX pixY
# returns: ulxNew ulyNew startCol startRow outRasCols outRasRows
# how to use it:
# use argument like "-projwin ulx uly lrx lry" in calling script
# use getRasterInfo function to get rasMinX rasMaxY pixX pixY in calling script
# use ReadAsArray in the calling script; give it startCol, startRow, outRasCols, outRasRows
# ... to perform your subset; you'll need in raster columns and rows from getRasterInfo too
# TODO: ADD SOMETHING THAT WILL SUBSET PARTIAL OVERLAP BETWEEN PROJWIN AND RASTER
def calcProjWin(ulx,uly,lrx,lry,rasMinX,rasMaxY,rasMaxX,rasMinY,rasPixX,rasPixY):
	[ulx, uly, lrx, lry,rasMinX,rasMaxY,rasMaxX,rasMinY,rasPixX,rasPixY] = [float(x) for x in [ulx, uly, lrx, lry,rasMinX,rasMaxY,rasMaxX,rasMinY,rasPixX,rasPixY]]
	# snap ulx, uly, lrx, lry to pixels of inRas; note that inRasPixY is -ve, so it is still subtracted just like inRasPixX
	if ulx > rasMinX:
		ulxNew = rasMinX + floor((ulx - rasMinX)/rasPixX) * rasPixX
	else:
		ulxNew = rasMinX
	if uly < rasMaxY:
		ulyNew = rasMaxY + floor((uly - rasMaxY)/rasPixY) * rasPixY
	else:
		ulyNew = rasMaxY
	if lrx < rasMaxX:
		lrxNew = rasMinX + ceil((lrx - rasMinX) / rasPixX) * rasPixX
		#lrxNew = rasMaxX + floor((rasMaxX - lrx)/rasPixX) * rasPixX
	else:
		lrxNew = rasMaxX
	if lry > rasMinY:
		lryNew = rasMaxY + ceil((lry - rasMaxY) / rasPixY) * rasPixY
		#lryNew = rasMinY + ceil((lry - rasMinY)/rasPixY) * rasPixY
	else:
		lryNew = rasMinY
	####### NEED TO UPDATE FROM HERE DOWN CONFIRM THAT IT WORKS WITH SUBSET.PY
	startCol = int( (ulxNew - rasMinX)/rasPixX )
	startRow = int( (ulyNew - rasMaxY)/rasPixY )
	outRasCols = int((lrxNew - ulxNew) / rasPixX)
	outRasRows = int((lryNew - ulyNew) / rasPixY)
	if outRasCols < 0:
		outRasCols = 0
	if outRasRows < 0:
		outRasRows = 0
	#print 'ulyNew: '+str(ulyNew)
	#print 'rasMaxY: '+str(rasMaxY)
	#print 'startCol '+str(startCol)+'startRow '+str(startRow)
	#outRasRows = outRasRows + 1    # add one pixel in Y direction
	#ulyNew = ulyNew - (1*rasPixY)  # adjust ULY to reflect addition in Y direction
	#outRasCols = outRasCols + 1    # add one pixel in X direction; no need to adjust ULX because pixel is added to the end
	
	#if startCol < 0:	startCol = 0
	#if startRow < 0:	startRow = 0
	#ulyNew = rasMaxY + ceil((uly - rasMaxY) / rasPixY) * rasPixY
	#lryNew = rasMaxY + ceil((lry - rasMaxY) / rasPixY) * rasPixY
	#startRow = int( (rasMaxY - ulyNew)/rasPixX )

	#ulxNew = rasMinX + floor((ulx - rasMinX) / rasPixX) * rasPixX
	#ulyNew = rasMaxY + floor((uly - rasMaxY) / rasPixY) * rasPixY
	#lrxNew = rasMinX + floor((lrx - rasMinX) / rasPixX) * rasPixX
	#lryNew = rasMaxY + floor((lry - rasMaxY) / rasPixY) * rasPixY
	return([ulxNew,ulyNew,lrxNew,lryNew,startCol,startRow,outRasCols,outRasRows])









