
from __future__ import division

import sys
import os

import lal
from lal import gpstime

from glue.ligolw import table as ligolw_table
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

from lalinference.bayestar import ligolw as ligolw_bayestar

file1 = sys.argv[1]
(path1, filename1) = os.path.split(file1)

file2 = sys.argv[2]
(path2, filename2) = os.path.split(file2)

xmldoc1 = ligolw_utils.load_filename(file1,
                                        contenthandler=ligolw_bayestar.LSCTablesContentHandler)

xmldoc2 = ligolw_utils.load_filename(file2,
                                        contenthandler=ligolw_bayestar.LSCTablesContentHandler)

sim_inspiral_table = ligolw_table.get_table(xmldoc1,
                                        lsctables.SimInspiralTable.tableName)

coinc_inspiral_table = ligolw_table.get_table(xmldoc2, 
                                              lsctables.CoincInspiralTable.tableName)

print "# GPS mass1 mass2 dist SNR RA dec inclination skymap"

for coinc in coinc_inspiral_table:

        for sim_inspiral in sim_inspiral_table:
            if abs(sim_inspiral.geocent_end_time-coinc.end_time) < 2:
                break
            else:
                sim_inspiral = None

        if not sim_inspiral:
            continue

        end_time = lal.LIGOTimeGPS(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
        (RA,dec) = sim_inspiral.ra_dec
        event_id = str(coinc.coinc_event_id)

        print '{gps};"{date}";{mass1};{mass2};{dist};{snr};{RA};{dec};{inclination};{skymap}'.format(gps=end_time,
                                                                                                     date=lal.gpstime.gps_to_utc(end_time).isoformat(' '),
                                                                                                     mass1=sim_inspiral.mass1,
                                                                                                     mass2=sim_inspiral.mass2,
                                                                                                     dist=sim_inspiral.distance,
                                                                                                     snr=coinc.snr,
                                                                                                     RA=RA,
                                                                                                     dec=dec,
                                                                                                     inclination=sim_inspiral.inclination,                     
                                                                                                     skymap="{}/{}.toa_phoa_snr.fits.gz".format(path1,event_id[-1]))
