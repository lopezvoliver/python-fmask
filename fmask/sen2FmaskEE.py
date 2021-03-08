
"""
Implementation of the cloud and shadow algorithms, specifically 
the detection of potential cloud and shadow pixels that form part of Fmask
for use in the Google Earth Engine (GEE) with Sentinel 2 (1C) data.

Currently, only the ee.Algorithms.FMask.fillMinima and ee.Algorithms.FMask.matchClouds
functions have been implemented in GEE 
https://developers.google.com/earth-engine/apidocs/ee-algorithms-fmask-fillminima
https://developers.google.com/earth-engine/apidocs/ee-algorithms-fmask-matchclouds

However, there is no straightforward method to obtain the required inputs
to the matchclouds algorithm, especially with Sentinel 2 data that does not
have a thermal band. 

Since the original python-fmask (fmask.py) functions are designed to work
with local rasters (using rios), there is no easy way to use these functions in GEE either.

This module is intended to solve this problem by offering functions that are
nearly identical to the fmask.py ones but with ee.Images as input

The cloud and shadow algorithms known collectively as Fmask are as published in 

Zhu, Z. and Woodcock, C.E. (2012). 
Object-based cloud and cloud shadow detection in Landsat imagery
Remote Sensing of Environment 118 (2012) 83-94. 
    
and
    
Zhu, Z., Wang, S. and Woodcock, C.E. (2015).
Improvement and expansion of the Fmask algorithm: cloud, cloud
shadow, and snow detection for Landsats 4-7, 8, and Sentinel 2 images
Remote Sensing of Environment 159 (2015) 269-277.
    
The original python-fmask was taken from Neil Flood's implementation by permission.

The notation and variable names are largely taken from the paper. Equation
numbers are also from the paper. 

Input is a top of atmosphere (TOA) GEE image (ee.Image - https://developers.google.com/earth-engine/apidocs/ee-image). 

The output is the potential cloud mask image, potential shadow mask image, and brightness temperature
image required for the ee.Algorithms.FMask.matchClouds algorithm.

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
from __future__ import print_function, division

import numpy
numpy.seterr(all='raise')
#from osgeo import gdal
#gdal.UseExceptions()
#from scipy.ndimage import uniform_filter, maximum_filter, label
#import scipy.stats

# our wrappers for bits of C that are installed with this package
#from . import fillminima
#from . import valueindexes
## configuration classes
from . import config
## exceptions
#from . import fmaskerrors
## so we can check if thermal all zeroes
#from . import zerocheck

# Bands in the saturation mask, if supplied
#SATURATION_BLUE = 0
#SATURATION_GREEN = 1
#SATURATION_RED = 2
    
#def potentialCloudFirstPass_ee(info, fmaskConfig, input_ee_img, outputs, otherargs):
def potentialCloudFirstPass_ee(fmaskConfig, input_ee_img, otherargs):

    """
    
    Calculate the first pass potential cloud layer (equation 6)
        
    """
    try:
        import ee 

    except Exception:
        raise ImportError(
            "You need to install the earth engine (ee) python API first."
        )

    ref = input_ee_img.divide(fmaskConfig.TOARefScaling)
    
    # Clamp off any reflectance <= 0
    #ref[ref<=0] = 0.00001  TODO: check if this is needed with S2 data coming from GEE

    # Extract the bands we need     
    # Reason for the "+1": fmaskConfig.bands[config.BAND_RED] for Sentinel 2 is 3 (e.g. B4 ). We need +1 to get "B"+"4"
    # The original code uses a number to select a band, with numbers starting from 0. e.g. 3 is the fourth element (0,1,2,3)
    # But here we select bands with a string such as "B4" so we need the actual number, not the position in an array
    # But for Swir1 and 2 it's ok as is -  e.g. B11 and B12
    blue = 'B{}'.format(fmaskConfig.bands[config.BAND_BLUE]+1)
    green = 'B{}'.format(fmaskConfig.bands[config.BAND_GREEN]+1)
    red = 'B{}'.format(fmaskConfig.bands[config.BAND_RED]+1)  
    nir = 'B{}'.format(fmaskConfig.bands[config.BAND_NIR]+1)
    swir1 = 'B{}'.format(fmaskConfig.bands[config.BAND_SWIR1])
    swir2 = 'B{}'.format(fmaskConfig.bands[config.BAND_SWIR2])
    #if hasattr(inputs, 'thermal'):   # TODO - for the moment, no thermal band in Sentinel 2
    #    THERM = otherargs.thermalInfo.thermalBand1040um
    #Note: otherargs.refBands is the same as fmaskConfig.bands 
    # config is already a module

    # Special mask needed only for resets in final pass
    # Not sure if any of this is needed:
    #refNullmask = (inputs.toaref[otherargs.bandsForRefNull] == otherargs.refNull).any(axis=0)

    #if hasattr(inputs, 'thermal'):
    #    thermNullmask = (inputs.thermal[THERM] == otherargs.thermalNull)
    #    nullmask = (refNullmask | thermNullmask)
    #    # Brightness temperature in degrees C
    #    bt = otherargs.thermalInfo.scaleThermalDNtoC(inputs.thermal)
    #else:
    #    thermNullmask = numpy.zeros_like(ref[0], dtype=numpy.bool)
    #    nullmask = refNullmask
    
    # Equation 1
    ndsi  = ref.normalizedDifference([green, swir1])
    ndvi = ref.normalizedDifference([nir,red]) 

    # In two parts, in case we have no thermal.
    basicTest = ref.select(swir2).gt(fmaskConfig.Eqn1Swir2Thresh).And(ndsi.lt(0.8)).And(ndvi.lt(0.8))
    #if hasattr(inputs, 'thermal'):  TODO - for the moment not needed for Sentinel 2
    #    basicTest = (basicTest & (bt < fmaskConfig.Eqn1ThermThresh))
    
    # Equation 2
    meanVis = (ref.select(blue).add(ref.select(green)).add(ref.select(red))).divide(3.0)
    
    whiteness = ((ref.select(blue).subtract(meanVis)).divide(meanVis)).abs()
    for n in [green, red]:
        whiteness = whiteness.add(((ref.select(n).subtract(meanVis)).divide(meanVis)).abs())
    
    whitenessTest = whiteness.lt(fmaskConfig.Eqn2WhitenessThresh)   

    # Equation 3 
    hazeTest = ((ref.select(blue).subtract(ref.select(red).multiply(0.5))).subtract(0.08)).gt(0)
    
    # Equation 4
    b45test = (ref.select(nir).divide(ref.select(swir1))).gt(0.75)
    
    # Equation 5
    waterTest = (ndvi.lt(0.01).And(ref.select(nir).lt(0.11))).Or(ndvi.lt(0.1).And(ref.select(nir).lt(0.05)))

    #waterTest[nullmask] = False #  TODO
    
    if config.BAND_CIRRUS in fmaskConfig.bands:
        # Zhu et al 2015, section 2.2.1. 
        cirrus = 'B{}'.format(fmaskConfig.bands[config.BAND_CIRRUS])
        cirrusBandTest = ref.select(cirrus).gt(fmaskConfig.cirrusBandTestThresh)
    
    # Equation 6. Potential cloud pixels (first pass)
    pcp = basicTest.And(whitenessTest).And(hazeTest).And(b45test)
    
    # If Sentinel-2, we can use the Frantz 2018 displacement test
    if (fmaskConfig.sensor == config.FMASK_SENTINEL2) and fmaskConfig.sen2displacementTest:
        cdi = ee.Algorithms.Sentinel2.CDI(input_ee_img)
        if otherargs.replace_and:
            selection = pcp.Or(cdi.lt(-0.5)) # I think this would work better with OR.. TODO add to config options
        else:
            selection = pcp.And(cdi.lt(-0.5)) # Default is to use this one 

        selection = selection.focal_min(otherargs.res)  # ee.Image.focal_min() is equivalent to scipy.ndimage.binary_erosion()
        # region grow within (cdi < -0.25)
        rg_mask = pcp.And(cdi.lt(-0.25))
        # for binary_dilation we use .focal_max
        # but it does not have the option iterations=0 as in scipy 
        # for now, we can use a few iterations (N_focal_max) 
        plus_k = ee.Kernel.plus(radius=otherargs.res, units='meters', normalize=True)
        selection_dilated = selection.focal_max(kernel=plus_k, iterations=otherargs.N_focal_max) 
        selection = selection.where(rg_mask, selection_dilated)
        pcp = pcp.updateMask(selection)
    

    # Include cirrusBandTest, from 2015 paper. Zhu et al. are not clear whether it is
    # supposed to be combined with previous tests using AND or OR, so I tried both
    # and picked what seemed best. 
    if config.BAND_CIRRUS in fmaskConfig.bands:
        pcp = (pcp.Or(cirrusBandTest))

    # TODO: check if we need this extra test 
    # This is an extra saturation test added by DERM, and is not part of the Fmask algorithm. 
    # However, some cloud centres are saturated, and thus fail the whiteness and haze tests
    #if hasattr(inputs, 'saturationMask'):
    #    saturatedVis = (inputs.saturationMask != 0).any(axis=0)
    #    veryBright = (meanVis > 0.45)
    #    saturatedAndBright = saturatedVis & veryBright
    #    pcp[saturatedAndBright] = True
    #    whiteness[saturatedAndBright] = 0
    
    #pcp[nullmask] = False
    
    # Equation 7
    clearSkyWater = waterTest.And(ref.select(swir2).lt(fmaskConfig.Eqn7Swir2Thresh))
    
    # Equation 12
    clearLand = pcp.Not().And(waterTest.Not())
    
    # Equation 15
    # Need to modify ndvi/ndsi by saturation......
    #if hasattr(inputs, 'saturationMask'): # TODO - check for the general case (no saturation mask in sentinel2Stacked)
    #    modNdvi = numpy.where((inputs.saturationMask[SATURATION_GREEN] != 0), 0, ndvi)
    #    modNdsi = numpy.where((inputs.saturationMask[SATURATION_RED] != 0), 0, ndsi)
    #else:
    #    modNdvi = ndvi
    #    modNdsi = ndsi

    modNdvi = ndvi
    modNdsi = ndsi

    # Maximum of three indices
    maxNdx = modNdvi.abs()
    maxNdx = maxNdx.max(modNdsi.abs())
    maxNdx = maxNdx.max(whiteness)
    variabilityProb = ee.Image([1]).subtract(maxNdx)

    #variabilityProb[nullmask] = 0
    #variabilityProbPcnt = numpy.round(variabilityProb * PROB_SCALE)
    #variabilityProbPcnt = variabilityProbPcnt.clip(BYTE_MIN, BYTE_MAX).astype(numpy.uint8)
    #variabilityProbPcnt - stored as 8 bits integer..not sure if needed
    # for now, the function will return the original variabilityProb
    
    # Equation 20
    # In two parts, in case we are missing thermal
    snowmask = (ndsi.gt(0.15).And(ref.select(nir).gt(fmaskConfig.Eqn20NirSnowThresh))).And(ref.select(green).gt(fmaskConfig.Eqn20GreenSnowThresh))
    #if hasattr(inputs, 'thermal'):  TODO - for now not needed for Sentinel 2
    #    snowmask = snowmask & (bt < fmaskConfig.Eqn20ThermThresh)
    #snowmask[nullmask] = False
    
    # Output the pcp and water test layers. 
    #outputs.pass1 = numpy.array([pcp, waterTest, clearLand, variabilityProbPcnt, 
    #    nullmask, snowmask, refNullmask, thermNullmask])
    return pcp, waterTest, clearLand, variabilityProb, snowmask

    # Accumulate histograms of temperature for land and water separately
    # TODO - not sure if this will be needed
    #if hasattr(inputs, 'thermal'):
    #    scaledBT = (bt + BT_OFFSET).clip(0, BT_HISTSIZE)
    #    otherargs.waterBT_hist = accumHist(otherargs.waterBT_hist, scaledBT[clearSkyWater])
    #    otherargs.clearLandBT_hist = accumHist(otherargs.clearLandBT_hist, scaledBT[clearLand])
    #scaledB4 = (ref[nir] * B4_SCALE).astype(numpy.uint8)
    #otherargs.clearLandB4_hist = accumHist(otherargs.clearLandB4_hist, scaledB4[clearLand])
    #otherargs.nonNullCount += numpy.count_nonzero(~nullmask)