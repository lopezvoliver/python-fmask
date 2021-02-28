# This file is part of 'python-fmask' - a cloud masking module
# Copyright (C) 2015  Neil Flood
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
Make a 4-layer image of satellite and sun angles, from the given tile metadata
file. 

The sun angles are exactly as provided in the XML. The satellite angles are more 
complicated, because they vary slightly with each band. I have done some limited 
testing and concluded that the variation is generally quite small, and so the 
layers written here are the per-pixel averages over all bands. The satellite
zenith this appears to vary only a couple of degrees across bands, while the 
variation for satellite azimuth is somewhat larger, depending on which bands. I
think that the 60m bands vary a bit more, but not sure. This is also complicated
by the variation across the scan, due to the different view angles of each detector
module. The pushbroom appears to be made up of 12 such modules, each looking
in slightly different directions, and each band within each module looking slightly
differently. Complicated...... sigh......

A bit of a description of how the instrument is structured can be found at
    https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-2-msi/msi-instrument

"""
from __future__ import print_function, division

import sys
import argparse

import numpy
from osgeo import gdal
gdal.UseExceptions()
from osgeo import osr

from rios import applier
from rios import calcstats
from rios import cuiprogress

from fmask import sen2meta

# This scale value will convert between DN and radians in output image file, 
#    radians = dn * SCALE_TO_RADIANS
SCALE_TO_RADIANS = 0.01


def getCmdargs():
    """
    Get commandline arguments
    """
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--infile", help="Input sentinel-2 tile metafile")
    p.add_argument("-o", "--outfile", help="Output angles image file")
    cmdargs = p.parse_args()
    if cmdargs.infile is None or cmdargs.outfile is None:
        p.print_help()
        sys.exit()
        
    return cmdargs


def mainRoutine():
    """
    Main routine
    """
    cmdargs = getCmdargs()
    
    makeAngles(cmdargs.infile, cmdargs.outfile)


def makeAngles(infile, outfile):
    """
    Callable main routine
    """
    info = sen2meta.Sen2TileMeta(filename=infile)
    
    ds = createOutfile(outfile, info)
    nullValDN = 1000
    
    # Get a sorted list of the Sentinel-2 band names. Note that sometimes this
    # is an incomplete list of band names, which appears to be due to a bug in 
    # earlier versions of ESA's processing software. I suspect it relates to 
    # Anomaly number 11 in the following page. 
    # https://sentinel.esa.int/web/sentinel/news/-/article/new-processing-baseline-for-sentinel-2-products
    bandNames = sorted(info.viewAzimuthDict.keys())
    
    # Mean over all bands
    satAzDeg = numpy.array([info.viewAzimuthDict[i] for i in bandNames])
    satAzDegMeanOverBands = satAzDeg.mean(axis=0)
    
    satZenDeg = numpy.array([info.viewZenithDict[i] for i in bandNames])
    satZenDegMeanOverBands = satZenDeg.mean(axis=0)

    sunAzDeg = info.sunAzimuthGrid
    
    sunZenDeg = info.sunZenithGrid
    
    stackDeg = numpy.array([satAzDegMeanOverBands, satZenDegMeanOverBands, sunAzDeg, sunZenDeg])
    stackRadians = numpy.radians(stackDeg)
    
    stackDN = numpy.round(stackRadians / SCALE_TO_RADIANS).astype(numpy.int16)
    nullmask = numpy.isnan(stackDeg)
    stackDN[nullmask] = nullValDN
    
    lnames = ['SatelliteAzimuth', 'SatelliteZenith', 'SunAzimuth', 'SunZenith']
    for i in range(ds.RasterCount):
        b = ds.GetRasterBand(i+1)
        b.WriteArray(stackDN[i])
        b.SetNoDataValue(nullValDN)
        b.SetDescription(lnames[i])
    calcstats.calcStats(ds, ignore=nullValDN, progress=cuiprogress.SilentProgress())
    del ds
        

def createOutfile(filename, info):
    """
    Create the empty output image file
    """
    drvr = gdal.GetDriverByName(applier.DEFAULTDRIVERNAME)
    (nrows, ncols) = info.anglesGridShape
    ds = drvr.Create(filename, ncols, nrows, 4, gdal.GDT_Int16, applier.DEFAULTCREATIONOPTIONS)
    gt = (info.anglesULXY[0], info.angleGridXres, 0, info.anglesULXY[1], 0.0, -info.angleGridYres)
    ds.SetGeoTransform(gt)
    
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(int(info.epsg))
    ds.SetProjection(sr.ExportToWkt())
    return ds

def makeAngles_ee(infile):
    """
    Callable function from python 
    to make the angles image as an ee.Image object
    infile is the meta (.xml) file and is read directly 
    from a URL (or file)
    This function is identical to the main function up until
    the writing to file part. Here we instead create the ee.Image to return 
    """
    import ee
    info = sen2meta.Sen2TileMeta(filename=infile)
    
    nullValDN = 1000
    
    # Get a sorted list of the Sentinel-2 band names. Note that sometimes this
    # is an incomplete list of band names, which appears to be due to a bug in 
    # earlier versions of ESA's processing software. I suspect it relates to 
    # Anomaly number 11 in the following page. 
    # https://sentinel.esa.int/web/sentinel/news/-/article/new-processing-baseline-for-sentinel-2-products
    bandNames = sorted(info.viewAzimuthDict.keys())
    
    # Mean over all bands
    satAzDeg = numpy.array([info.viewAzimuthDict[i] for i in bandNames])
    satAzDegMeanOverBands = satAzDeg.mean(axis=0)
    
    satZenDeg = numpy.array([info.viewZenithDict[i] for i in bandNames])
    satZenDegMeanOverBands = satZenDeg.mean(axis=0)

    sunAzDeg = info.sunAzimuthGrid
    
    sunZenDeg = info.sunZenithGrid
    
    stackDeg = numpy.array([satAzDegMeanOverBands, satZenDegMeanOverBands, sunAzDeg, sunZenDeg])
    stackRadians = numpy.radians(stackDeg)
    
    stackDN = numpy.round(stackRadians / SCALE_TO_RADIANS).astype(numpy.int16)
    nullmask = numpy.isnan(stackDeg)
    stackDN[nullmask] = nullValDN
    
    lnames = ['SatelliteAzimuth', 'SatelliteZenith', 'SunAzimuth', 'SunZenith']

    # Write data to image. 
    # see Q here: https://groups.google.com/g/google-earth-engine-developers/c/JKO3YWKMaAw/m/_nqBju2XAQAJ
    # and A here: https://code.earthengine.google.com/9a04e6f575444bc84dd4fcddf5fdc2ea
    
    # We start with the first band (satellite azimuth):
    data_to_write = ee.Array([
        stackDN[0].tolist()
    ]).reshape(list(info.anglesGridShape)).transpose()
    # Geotransform correct order for ee (not the same as in gdal!):
    #[xScale, xShearing, xTranslation, yShearing, yScale, yTranslation]
    gt = [float(info.angleGridXres), float(0), float(info.anglesULXY[0]), float(0), float(-info.angleGridYres), float(info.anglesULXY[1])]
    projection = ee.Projection('EPSG:{}'.format(info.epsg), gt)
    # Construct an image of (x,y) pixel coordinates in our projection.
    coords = ee.Image.pixelCoordinates(projection).floor().int32()
    # Clamp the coordinates to the valid range of indices into the data array.
    dimx = data_to_write.length().getInfo()[0]
    dimy = data_to_write.length().getInfo()[1]
    x = coords.select('x')
    y = coords.select('y')
    coords = coords.updateMask(x.gte(0).And(y.gte(0)).And(x.lt(dimx)).And(y.lt(dimy)))
    # Now index into the array at the computed coordinates.
    satAzRadMeanOverBands_ee = ee.Image(data_to_write).arrayGet(coords)

    # Repeat for the other bands:
    data_to_write = ee.Array([
        stackDN[1].tolist()
    ]).reshape(list(info.anglesGridShape)).transpose()
    satZenRadMeanOverBands_ee = ee.Image(data_to_write).arrayGet(coords)
    
    data_to_write = ee.Array([
        stackDN[2].tolist()
    ]).reshape(list(info.anglesGridShape)).transpose()
    sunAzRad_ee = ee.Image(data_to_write).arrayGet(coords)
    
    data_to_write = ee.Array([
        stackDN[3].tolist()
    ]).reshape(list(info.anglesGridShape)).transpose()
    sunZenRad_ee = ee.Image(data_to_write).arrayGet(coords)
    
    # Create the image by appending the bands:
    angles_img = satAzRadMeanOverBands_ee.addBands(satZenRadMeanOverBands_ee).addBands(sunAzRad_ee ).addBands(sunZenRad_ee)
    angles_img = angles_img.select(['constant','constant_1','constant_2', 'constant_3'],lnames)
    
    # update mask to not show the null values (1000)?
    angles_img = angles_img.updateMask(angles_img.neq(nullValDN))

    return(angles_img)
    
    
    