from __future__ import division

import sys
import math
import os
from shutil import copy

from glue.ligolw import ligolw
from glue.ligolw import ilwd
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils

LITTLEHOPE_HOME = os.getcwd()
LITTLEHOPE_OPTS = '--detector H1 --detector L1 --min-triggers 2 --snr-threshold 4.0 --reference-psd {home}/psds_2016-17.xml --template-bank {home}/templates.xml --waveform "TaylorF2threePointFivePN"'.format(home=LITTLEHOPE_HOME)
CMD_LITTLEHOPE = '{home}/my_bayestar_littlehope  {opts} {simdir}/{mdc_file} -o coinc.xml\n'

LOCALCOINCS_OPTS = '--waveform "TaylorF2threePointFivePN" --f-low 30'
CMD_LOCALCOINCS = 'bayestar_localize_coincs {opts} coinc.xml\n'
CMD_PLOT = 'for f in *.fits.gz; do bayestar_plot_allsky $f --contour 90 --radec 0.0 0.0 -o ${f%.*}.png; done\n'
CMD_JOBSUB = '/opt/sge/bin/lx24-amd64/qsub -N {} -o {} -e {} {} {}'

if __name__ == "__main__":

    # parse input args
    if len(sys.argv) <= 1:
        sys.exit()

    indir = sys.argv[1]
    
    # list contents of input dir
    if os.path.isdir(indir):
        files = os.listdir(indir)
    else:
        sys.exit()

    if len(sys.argv) > 2:
        jobsub_opts = sys.argv[2]
    else:
        jobsub_opts = ""

    # loop on coincX.xml
    submission_cmds = []
    select_mdc_xml = lambda f: f.lower().startswith('mdc') and f.lower().endswith('.xml')
    for file in filter(select_mdc_xml, files):
        
        # create job dirs
        simdir,_ = os.path.splitext(file)
        if not os.path.exists(simdir):
            os.makedirs(simdir)

        copy(indir + '/' + file, simdir + '/')
        
        fullpath_simdir = os.path.abspath(simdir)

        # create job scripts
        batch_filename = "{}/batch.sh".format(simdir)
        with open(batch_filename, "w") as batch_script:
            batch_script.write("#!/usr/bin/env bash\n")
            batch_script.write("cd {simdir}\n".format(simdir=fullpath_simdir))
            batch_script.write(CMD_LITTLEHOPE.format(mdc_file=file,
                                                     simdir=fullpath_simdir,
                                                     home=LITTLEHOPE_HOME,
                                                     opts=LITTLEHOPE_OPTS))
            batch_script.write(CMD_LOCALCOINCS.format(simdir=fullpath_simdir,
                                                      opts=LOCALCOINCS_OPTS))
            batch_script.write(CMD_PLOT)
            batch_script.write('echo "hello word!"')

        os.chmod(batch_filename, 0744)
        
        # create script with submission list
        submission_cmds.append(CMD_JOBSUB.format(simdir, 
                                                 "{}/{}.out".format(fullpath_simdir,simdir), 
                                                 "{}/{}.err".format(fullpath_simdir,simdir), 
                                                 jobsub_opts,
                                                 batch_filename))

    outfile = open('jobsubmission.sh', 'w')
    print>>outfile, "#!/usr/bin/env bash"
    for line in submission_cmds:
        print>>outfile, line

    os.chmod('jobsubmission.sh', 0744)
