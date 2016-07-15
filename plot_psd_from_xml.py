from __future__ import division

import sys
import math
import numpy

import matplotlib
import matplotlib.pyplot as plot

import lal

from lal.lal import LIGOTimeGPS
import lal.series

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import table as ligolw_table
from glue.ligolw import lsctables

import timing
from lalinference.bayestar import filter

font = {'family' : 'helvetica',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
        pass

# from the Prospect paper
#
# data = numpy.array([
#     [10, 1.0e-43,  numpy.nan],
#     [20, 8.1e-45,  numpy.nan],
#     [30, 1.6e-45,   3.6e-43],
#     [40, numpy.nan,   4e-44],
#     [60,  1.2e-46,  2.5e-45],
#     [100, 4.9e-47,  4.4e-46],
#     [200, 4.9e-47,  4.4e-46],
#     [500, 8.1e-47,  1.0e-45],
#     [1000, numpy.nan, 3.6e-45],
#     [1200, 4e-46,  numpy.nan],
#     [2000, 9e-46,   1.4e-44]])

data_adv = numpy.loadtxt("Adv_Virgo_20Mpc_PSD.txt",dtype={'names':('f','psd'), 'formats':('f8','f8')})
data_aligo = numpy.loadtxt("aLIGO_80Mpc_PSD.txt",dtype={'names':('f','psd'), 'formats':('f8','f8')})

filename = "psds_2016-17.xml"
#filename = "psd.xml"
duration = 32 # sec
sample_rate = 1024 # Hz

xmldoc = ligolw_utils.load_filename(
    filename, contenthandler=lal.series.PSDContentHandler)
input_psds = lal.series.read_psd_xmldoc(xmldoc)
psd_models = dict(
    (key, timing.InterpolatedPSD(filter.abscissa(psd), psd.data.data))
    for key, psd in input_psds.iteritems() if psd is not None)

for detector, model in psd_models.iteritems():
    psd = lal.CreateREAL8FrequencySeries(None, lal.LIGOTimeGPS(0), 0, 1./duration, filter.unitInverseHertz, int((duration * sample_rate) / 2 + 1))
    psd.data.data = model(filter.abscissa(psd))
    plot.loglog(filter.abscissa(psd), psd.data.data, label=detector)

# plot.loglog(data[:,0],data[:,1])
# plot.loglog(data[:,0],data[:,2])

plot.loglog(data_adv['f'],data_adv['psd'], label="V1 from Prospects paper")
plot.loglog(data_aligo['f'],data_aligo['psd'], label="H1 and L1 from Prospects paper")

plot.legend(loc=1)
plot.grid(True)
plot.xlabel('Frequency Hz')
plot.ylabel('Strain 1/Hz')
plot.show()
