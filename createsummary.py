
from __future__ import division

import sys
import lal
from lal import gpstime

from glue.ligolw import table as ligolw_table
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

from lalinference.bayestar import ligolw as ligolw_bayestar

xmldoc1 = ligolw_utils.load_filename(sys.argv[1],
                                        contenthandler=ligolw_bayestar.LSCTablesContentHandler)
sim_inspiral_table = ligolw_table.get_table(xmldoc1,
                                        lsctables.SimInspiralTable.tableName)

xmldoc2 = ligolw_utils.load_filename(sys.argv[2],
                                         contenthandler=ligolw_bayestar.LSCTablesContentHandler)

print "# GPS mass1 mass2 dist SNR RA dec"
    
for coinc, sngl_inspirals in ligolw_bayestar.coinc_and_sngl_inspirals_for_xmldoc(xmldoc2):

    sngl_inspiral = sngl_inspirals[0]

    for sim_inspiral in sim_inspiral_table:
        if sim_inspiral.geocent_end_time == sngl_inspiral.end_time:
            break

    end_time = lal.LIGOTimeGPS(sim_inspiral.geocent_end_time, sim_inspiral.geocent_end_time_ns)
    (RA,dec) = sim_inspiral.ra_dec

    print "{gps};{date},{mass1};{mass2};{dist};{snr};{RA};{dec}".format(gps=end_time,
                                                                        date=lal.gpstime.gps_to_utc(end_time).isoformat(' '),
                                                                        mass1=sim_inspiral.mass1,
                                                                        mass2=sim_inspiral.mass2,
                                                                        dist=sim_inspiral.distance,
                                                                        #dist=sngl_inspiral.eff_distance,
                                                                        snr=sngl_inspiral.snr,
                                                                        RA=RA,
                                                                        dec=dec)
