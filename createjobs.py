from __future__ import division

import sys
import math
import os
from shutil import copy

from glue.ligolw import ligolw
from glue.ligolw import ilwd
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils

LITTLEHOPE_OPTS = '--detector H1 --detector L1 --min-triggers 2 --snr-threshold 4.0 --reference-psd psds_2016-17.xml --template-bank templates.xml --waveform "TaylorF2threePointFivePN"'
CMD_LITTLEHOPE = './my_bayestar_littlehope  {} {{}} -o coinc.xml'.format(LITTLEHOPE_OPTS)

LOCALCOINCS_OPTS = '--waveform "TaylorF2threePointFivePN" --f-low 30'
CMD_LOCALCOINCS = 'bayestar_localize_coincs {} coinc.xml'.format(LOCALCOINCS_OPTS)
CMD_JOBSUB = '/opt/sge/bin/lx24-amd64/qsub -N {} -o {} {} {}'

if __name__ == "__main__":

    # parse input args
    if len(sys.argv) > 0:
        indir = sys.argv[1]
    
    # list contents of input dir
    if os.path.isdir(indir):
        files = os.listdir(indir)
    else:
        sys.exit()

    # loop on coincX.xml
    submission_cmds = []
    select_mdc_xml = lambda f: f.lower().startswith('mdc') and f.lower().endswith('.xml')
    for file in filter(select_mdc_xml, files):
        
        # create job dirs
        simdir,_ = os.path.splitext(file)
        if not os.path.exists(simdir):
            os.makedirs(simdir)

        copy(indir + '/' + file, simdir + '/')
        
        # create job scripts
        batch_filename = "{}/batch.sh".format(simdir)
        with open(batch_filename, "w") as batch_script:
            batch_script.write("#!/usr/bin/env bash")
            batch_script.write("cd {}".format(simdir))
            batch_script.write(CMD_LITTLEHOPE.format(file))
            batch_script.write(CMD_LOCALCOINCS.format(file))
            #batch_script.write(CMD_PLOT.format(file))

        os.chmod(batch_filename, 0744)
        
        # create script with submission list
        submission_cmds.append(CMD_JOBSUB.format(simdir, "{}/{}.txt".format(simdir,simdir), "", batch_filename))

    outfile = open('jobsubmission.sh', 'w')
    print>>outfile, "#!/usr/bin/env bash"
    for line in submission_cmds:
        print>>outfile, line
