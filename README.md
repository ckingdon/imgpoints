> I wrote this circa 2010. Someday I'll re-write it, apply proper programming principles, and make it more efficient. Until then, here it is, warts and all!

# imgpoints.py
Command-line tool to quickly extract data at xy points from raster datasets.

```CONSOLE
usage: imgpoints.py [-h] [-p [POINT_INPUT]] [-r [RASTER_INPUT]]
                    [-c [CSV_OUTPUT]] [-x [XCOORD_FIELD]] [-y [YCOORD_FIELD]]
                    [-f [FID_FIELD]] [-e [EPSG_CODE]] [-k] [-n [NODATA_VALUE]]
                    [-i] [-t [TRANSFORMATION_FILE]] [-4]

extract pixel values based on point locatons

optional arguments:
  -h, --help            show this help message and exit
  -p [POINT_INPUT], --point-input [POINT_INPUT]
                        input point file; either .shp or .csv
  -r [RASTER_INPUT], --raster-input [RASTER_INPUT]
                        input raster file; anything supported by GDAL
  -c [CSV_OUTPUT], --csv-output [CSV_OUTPUT]
                        csv file created to hold extracted pixel values
  -x [XCOORD_FIELD], --xcoord-field [XCOORD_FIELD]
                        optional: point file's X coordinate field; only used
                        if POINT_INPUT is a .csv
  -y [YCOORD_FIELD], --ycoord-field [YCOORD_FIELD]
                        optional: point file's Y coordinate field; only used
                        if POINT_INPUT is a .csv
  -f [FID_FIELD], --fid-field [FID_FIELD]
                        optional: point file field containing UNIQUE
                        identifier; if not specified, a field named ZID will
                        be created in CSV_OUTPUT
  -e [EPSG_CODE], --epsg-code [EPSG_CODE]
                        optional: input EPSG code for point file; .prj
                        associated with .shp takes precedence
  -k, --keep-all        NOT WORKING YET!! optional: default=yes; keep all
                        point even if they don't intersect raster
  -n [NODATA_VALUE], --nodata-value [NODATA_VALUE]
                        NOT WORKING YET!! optional: value assigned to points
                        that don't intersect RASTER_FILE
  -i, --ignore-rotation
                        optional: ignore rotation in ENVI header
  -t [TRANSFORMATION_FILE], --transformation-file [TRANSFORMATION_FILE]
                        optional: output from wktproj2coeffile.py; modified
                        ENVI polynomial coefficient transformation file
  -4, --four-nearest    optional:collect values for 4 nearest pixels to point

depends on the following python modules: 
osgeo (gdal, ogr, osr, gdalconst), funcOSGEOproc.py (home-rolled)             
os, os.path, sys, argparse, csv, re
```
