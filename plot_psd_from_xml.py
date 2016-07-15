from __future__ import division

import sys
import math

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

#filename = "psds_2016-17.xml"
filename = "psd.xml"
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

plot.legend(loc=1)
plot.grid(True)
plot.xlabel('Frequency Hz')
plot.ylabel('Strain 1/Hz')
plot.show()
