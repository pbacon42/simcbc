
from __future__ import division

import sys
import math

import numpy
import numpy.random as rand

import lal
from lal import CachedDetectors
from lal import LALDetectorIndexLHODIFF,LALDetectorIndexLLODIFF,LALDetectorIndexVIRGODIFF

from lal.lal import LIGOTimeGPS
from lal.lal import MSUN_SI as LAL_MSUN_SI
from lal.lal import PC_SI as LAL_PC_SI
from lal.lal import DimensionlessUnit
from lal.lal import CreateREAL8FrequencySeries, CreateCOMPLEX16FrequencySeries

import lal.series

import lalsimulation
from lalsimulation.lalsimulation import SimInspiralTD, SimInspiralFD
from lalsimulation.lalsimulation import SimNoisePSD
from lalsimulation.lalsimulation import SimNoisePSDaLIGOZeroDetHighPower, SimNoisePSDVirgo
from lalsimulation.lalsimulation import SimDetectorStrainREAL8TimeSeries
from lalsimulation.lalsimulation import SimInspiralCreateWaveformFlags
from lalsimulation.lalsimulation import GetApproximantFromString

from pylal import antenna

from glue.ligolw import ligolw
from glue.ligolw import ilwd
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils

import timing
from lalinference.bayestar import filter

DETECTOR_SITES = {
    'H1': LALDetectorIndexLHODIFF,
    'L1': LALDetectorIndexLLODIFF,
    'V1': LALDetectorIndexVIRGODIFF
    }

DETECTOR_NOISE_MODELS = {
    'H1': SimNoisePSDaLIGOZeroDetHighPower,
    'L1': SimNoisePSDaLIGOZeroDetHighPower,
    'V1': SimNoisePSDVirgo
    }

DETECTOR_PSD_FILES = {
    'H1': "aLIGO_80Mpc_PSD_extended.txt",
    'L1': "aLIGO_80Mpc_PSD_extended.txt",
    'V1': "Adv_Virgo_20Mpc_PSD.txt"
    }

ZERO_SPIN = {'x': 0., 'y': 0., 'z': 0.}

# map order integer to a string that can be parsed by lalsimulation
PN_ORDERS = {
    'default'          : -1,
    'zeroPN'           : 0,
    'onePN'            : 2,
    'onePointFivePN'   : 3,
    'twoPN'            : 4,
    'twoPointFivePN'   : 5,
    'threePN'          : 6,
    'threePointFivePN' : 7,
    'pseudoFourPN'     : 8,
    }

START_O2 = 1156723217 # Thu Sep 01 00:00:00 GMT 2016
STOP_O2 = 1172361617   # Wed Mar 01 00:00:00 GMT 2017

class CompactBinary(object):
    """
    A CompactBinary object characterises a binary formed by two compact objects.
    """

    def __init__(self, mass1, mass2, distance, redshift, spin1, spin2, lambda1, lambda2, iota):
        """
        mass1, mass2 -- masses of the binary components in solar masses
        distance -- distance of the binary in Mpc
        redshift -- redshift of the binary. If zero, cosmology is ignored.
        spin1, spin2 -- spin vectors of binary components
        lambda1, lambda2 -- dimensionless tidal parameter of binary components
        iota -- inclination angle with respect to the line of sight in degrees
        """
        self.mass1 = mass1
        self.mass2 = mass2
        self.mchirp, self.eta = _mass1_mass2_to_mchirp_eta(mass1, mass2)
        self.distance = distance
        self.z = redshift        
        self.spin1 = spin1
        self.spin2 = spin2
        self.lambda1 = lambda1 
        self.lambda2 = lambda2
        self.iota = iota

class CBCTemplate(object):
    """
    A CBCTemplate object characterises the gravitational
    wave (GW) chirp signal associated to the coalescence of two
    inspiralling compact objects.
    """

    def __init__(self, approximant, amplitude0, phase0, sampling_rate, segment_duration, freq_min, freq_max, freq_ref, phi_ref, nonGRparams):
        """
        approximant -- model approximant
        amplitude0  -- amplitude pN order: -1 means include all
        phase0      -- phase pN order: -1 means include all
        sampling_rate    -- sampling rate in Hz
        segment_duration -- segment duration in sec
        freq_min -- start frequency in Hz
        freq_max -- end frequency in Hz
        freq_ref -- reference frequency for precessing spins in Hz
        phi_ref  -- final phase in degrees
        nonGRparams -- non GR parameters
        """
        
        self.approximant = GetApproximantFromString(approximant)
        self.sampling_rate = sampling_rate # Hz
        self.segment_duration = segment_duration # sec
        self.amplitude0 = amplitude0
        self.phase0 = phase0
        self.freq_min = freq_min # Hz, start frequency
        self.freq_max = freq_max # Hz, end frequency
        self.freq_ref = freq_ref # Hz, reference frequency for precessing spins
        self.phi_ref  = phi_ref  # final phase in degrees
        self.nonGRparams = nonGRparams # non GR parameters
        self.waveform_flags = SimInspiralCreateWaveformFlags()
        
    def freq_template(self, binary):
        """
        Computes the frequency-domain template model of the gravitational wave for a given compact binary.
        """
        
        frequency_resolution = 1.0 / self.segment_duration

        return SimInspiralFD(math.radians(self.phi_ref), frequency_resolution,
                                 binary.mass1 * LAL_MSUN_SI, binary.mass2 * LAL_MSUN_SI,
                                 binary.spin1['x'], binary.spin1['y'], binary.spin1['z'],
                                 binary.spin2['x'], binary.spin2['y'], binary.spin2['z'],
                                 self.freq_min, self.freq_max, self.freq_ref,
                                 binary.distance * 1.0e6 * LAL_PC_SI, binary.z, math.radians(binary.iota), binary.lambda1, binary.lambda2,
                                 self.waveform_flags, self.nonGRparams, self.amplitude0, self.phase0, self.approximant)

    def time_template(self, binary):
        """
        Compute time-domain template model of the gravitational wave for a given compact binary.
        """
    
        return SimInspiralTD(math.radians(self.phi_ref), 1.0 / self.sampling_rate,
                                 binary.mass1 * LAL_MSUN_SI, binary.mass2 * LAL_MSUN_SI,
                                 binary.spin1['x'], binary.spin1['y'], binary.spin1['z'],
                                 binary.spin2['x'], binary.spin2['y'], binary.spin2['z'],
                                 self.freq_min, self.freq_ref,
                                 binary.distance * 1.0e6 * LAL_PC_SI, binary.z,
                                 math.radians(binary.iota), binary.lambda1, binary.lambda2,
                                 self.waveform_flags, self.nonGRparams, self.amplitude0, self.phase0, self.approximant)
                                     
    # hstrain = SimDetectorStrainREAL8TimeSeries(hplus, hcross, ra, dec, psi, det)
    # hstrain.epoch += time_at_coalescence # set end time to time_at_coalescence
    # times = time_at_coalescence + sampling_period * numpy.arange(hstrain.data.length)  
    # signal = hstrain.data.data

class Detector(object):
    """
    A Detector object characterises a gravitational wave (GW) interferometric detector
    """

    def __init__(self, detector):
        """
        detector  -- label string of the detector
        descriptor -- LAL descriptor
        location -- geographic location of the detector
        response -- response matrix

        """
        self.name = detector
        self.descriptor =  CachedDetectors[DETECTOR_SITES[detector]]
        self.location = lalsimulation.DetectorPrefixToLALDetector(detector).location
        self.response = lalsimulation.DetectorPrefixToLALDetector(detector).response
        
    def antenna_pattern(self, time_at_coalescence, RA, dec, iota, psi):
        """ Compute antenna response
        """
        fplus,fcross,_,_ = antenna.response(time_at_coalescence,
                                    RA, dec, iota, psi, 'degree', self.name)
        return fplus, fcross
        
    def project_strain(self, hplus, hcross, time_at_coalescence, RA, dec, iota, psi):
        """ Project hplus and hcross onto the detector assuming a given
        position and polarization of the source.
        """

        assert hplus.data.length == hcross.data.length
        assert hplus.deltaF == hcross.deltaF
        assert hplus.f0 == hcross.f0

        freq_resolution = hplus.deltaF
        length = hplus.data.length

        fplus, fcross = self.antenna_pattern(time_at_coalescence, RA, dec, iota, psi)
    
        hstrain = CreateCOMPLEX16FrequencySeries("strain", 0.0, 0.0,
                                    freq_resolution,DimensionlessUnit,length);
        hstrain.data.data = fplus * hplus.data.data + fcross * hcross.data.data

        return hstrain

    def psd(self, freq_min, freq_resolution, length):
            """ Compute PSD from noise model
            """
            
            power_spec_density = CreateREAL8FrequencySeries("spectrum", 0.0, freq_min,
                                        freq_resolution, DimensionlessUnit, length);
            # SimNoisePSD(power_spec_density, freq_min, DETECTOR_NOISES[self.name])
            data = numpy.loadtxt(DETECTOR_PSD_FILES[self.name],dtype={'names':('freq','psd'), 'formats':('f8','f8')})
            model = timing.InterpolatedPSD(data['freq'], data['psd'])
            power_spec_density.data.data = model(filter.abscissa(power_spec_density))
            
            return power_spec_density

    def effective_distance(self, distance, time_at_coalescence, RA, dec, iota, psi):
        """ Returns the effective distance
        """
        
        fplus, fcross = self.antenna_pattern(time_at_coalescence, RA, dec, iota, psi)
        return distance / math.sqrt(fplus**2 + fcross**2)

    def time_delay_from_earth_center(self, RA, dec, time_gps):
        """ Returns the time delay from the earth center
        """
        return lal.TimeDelayFromEarthCenter(self.location,
                      float(RA), float(dec), float(time_gps))
 
H1 = Detector("H1")
L1 = Detector("L1")
    
def signal_to_noise(hstrain, psd, freq_min, freq_max):
    """
    Compute the signal-to-noise ratio of signal hstrain in detector noise
    in a specified frequency band.
    """

    # interpolate PSD on the same frequency axis
    hstrain_freqs = hstrain.f0 + hstrain.deltaF * numpy.arange(hstrain.data.length)
    psd_freqs = psd.f0 + psd.deltaF * numpy.arange(psd.data.length)
    psd_interp = numpy.interp(hstrain_freqs, psd_freqs, psd.data.data)

    selected_idx = (hstrain_freqs >= freq_min) & (hstrain_freqs <= freq_max)
        
    return math.sqrt(4.0 * hstrain.deltaF * numpy.sum(numpy.abs(hstrain.data.data[selected_idx])**2/psd_interp[selected_idx]))

def _empty_row(obj):
    """Create an empty sim_inspiral or sngl_inspiral row where the columns have
    default values of 0.0 for a float, 0 for an int, '' for a string. The ilwd
    columns have a default where the index is 0.
    """

    # check if sim_inspiral or sngl_inspiral
    if obj == lsctables.SimInspiral:
        row = lsctables.SimInspiral()
        cols = lsctables.SimInspiralTable.validcolumns
    else:
        row = lsctables.SnglInspiral()
        cols = lsctables.SnglInspiralTable.validcolumns

    # populate columns with default values
    for entry in cols.keys():
        if cols[entry] in ['real_4','real_8']:
            setattr(row,entry,0.)
        elif cols[entry] == 'int_4s':
            setattr(row,entry,0)
        elif cols[entry] == 'lstring':
            setattr(row,entry,'')
        elif entry == 'process_id':
            row.process_id = ilwd.ilwdchar("sim_inspiral:process_id:0")
        elif entry == 'simulation_id':
            row.simulation_id = ilwd.ilwdchar("sim_inspiral:simulation_id:0")
        elif entry == 'event_id':
            row.event_id = ilwd.ilwdchar("sngl_inspiral:event_id:0")
        else:
            raise ValueError("Column %s not recognized." %(entry) )

    return row

def _mass1_mass2_to_mchirp_eta(mass1, mass2):
    """ Convert mass1, mass2 into mchirp, eta params
    """
    m_total = mass1 + mass2
    eta = (mass1 * mass2) / (m_total * m_total)
    m_chirp = m_total * eta**(3./5.)
    return m_chirp,eta

if __name__ == "__main__":

    # parameters

    detectors = [H1, L1]
    
    time_from_start = 0  # s
    stride = 3600   # s
    jitter = 600    # s
    threshold = 4 # SNR selection threshold

    approximant = "TaylorT4treePN"
    amplitude_order = 0
    phase_order = -1
    sampling_rate = 1024 # Hz
    segment_duration = 64 # s
    freq_min = 10 # Hz
    freq_max = sampling_rate/2.0

    # create new sim and sngl tables
    
    class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
        pass
    
    lsctables.use_in(LIGOLWContentHandler)
    sim_table = lsctables.New(lsctables.SimInspiralTable)
    sngl_table = lsctables.New(lsctables.SnglInspiralTable,
                            columns=lsctables.SnglInspiralTable.validcolumns)

    counter = 0

    # loop through input file
    with open(sys.argv[1]) as infile:
      
        for line in infile:

            # ignore comment lines
            if line.startswith('##'):
                continue

            # read and parse one line
            try:
                distance, redshift, mass1, mass2, RA, dec = [float(x) for x in line.strip().split(" ")]
            except ValueError:
                continue

            # initialize lists of computed vars
            SNRs = []
            eff_distances = []
            end_times_at_detector = []
            PSDs = {}

            # generate end_time and remaining angle randomly
            geocent_end_time = START_O2 + time_from_start + rand.uniform(-jitter/2,jitter/2)
            iota = 360.0 * rand.random()  # all angles are in degrees
            phi_ref = 360.0 * rand.random()
            psi = 360.0 * rand.random()

            binary = CompactBinary(mass1, mass2, distance, redshift,
                                     ZERO_SPIN, ZERO_SPIN, 0.0, 0.0, iota)
            model = CBCTemplate(approximant, amplitude_order, phase_order, sampling_rate,
                                    segment_duration, freq_min, freq_max, 0.0, phi_ref, None)

            # compute strain model
            hplus, hcross = model.freq_template(binary)

            for detector in detectors:
                
                # project strain onto detector
                hstrain = detector.project_strain(hplus, hcross, geocent_end_time, RA, dec, iota, psi)

                # compute PSD -- XXX in principle, this should be computed and store only once XXX
                PSDs[detector.name] = detector.psd(model.freq_min,hstrain.deltaF,hstrain.data.length)
            
                # compute SNR
                SNRs.append(signal_to_noise(hstrain,
                                    PSDs[detector.name],
                                    model.freq_min, model.freq_max))

                # compute effective distance
                eff_distances.append(detector.effective_distance(distance, geocent_end_time,
                                                       RA, dec, iota, psi))

                # compute end time at detector
                time_delay = detector.time_delay_from_earth_center(RA, dec, geocent_end_time)
                end_times_at_detector.append(geocent_end_time + time_delay)


            # select injection if sufficient SNR at one of the detectors
            if all(snr < threshold for snr in SNRs):
                continue

            print "{} -- m1={} Msun m2={} Msun d={} Mpc -- SNR = {}".format(counter,
                                            binary.mass1, binary.mass2, binary.distance, SNRs)
        
            # create sim entry
            sim = _empty_row(lsctables.SimInspiral)
            
            sim.f_lower = model.freq_min
            sim.geocent_end_time = int(geocent_end_time)
            sim.geocent_end_time_ns = int(geocent_end_time % 1 * 1e9)
            sim.inclination = math.radians(iota)
            sim.latitude = math.radians(dec)
            sim.longitude = math.radians(RA)
            sim.mass1 = binary.mass1
            sim.mass2 = binary.mass2
            sim.mchirp = binary.mchirp
            sim.eta = binary.eta
            sim.polarization = math.radians(psi)
            sim.taper = 'TAPER_STARTEND'
            sim.distance = binary.distance
            sim.numrel_data = ""
            sim.spin1x = 0.0
            sim.spin1y = 0.0
            sim.spin1z = 0.0
            sim.spin2x = 0.0
            sim.spin2y = 0.0
            sim.spin2z = 0.0
            for det, end_time, eff_dist  in zip(detectors, end_times_at_detector, eff_distances):
                setattr(sim, det.name[0].lower()+'_end_time', int(end_time))
                setattr(sim, det.name[0].lower()+'_end_time_ns', int(end_time % 1 * 1e9))
                setattr(sim, 'eff_dist_'+det.name[0].lower(), eff_dist)

            # construct waveform string that can be parsed by lalsimulation
            #waveform_string = approximant
            # phase_order = lalsimulation.GetOrderFromString(s)
            #if not PN_ORDERS[phase_order] == -1:
            #    waveform_string += opts.order
            sim.waveform = approximant

            # append sim entry
            sim_table.append(sim)

            # create sngl entry
            sngl = _empty_row(lsctables.SnglInspiral)
            sngl.mass1 = binary.mass1
            sngl.mass2 = binary.mass2
            sngl.mchirp = binary.mchirp
            sngl.eta = binary.eta
            sngl.mtotal = binary.mass1 + binary.mass2
            sngl.spin1x = 0.0
            sngl.spin1y = 0.0
            sngl.spin1z = 0.0
            sngl.spin2x = 0.0
            sngl.spin2y = 0.0
            sngl.spin2z = 0.0
            sngl.end_time = int(geocent_end_time)
            sngl.end_time_ns = int(geocent_end_time % 1 * 1e9)

            # append sngl entry
            sngl_table.append(sngl)

            # increment time for next injection
            time_from_start += stride
            counter += 1

        print "{} mergers selected".format(counter)
        
        # build and write injection XML document
        xmldoc_mdc = ligolw.Document()
        xmldoc_mdc.appendChild(ligolw.LIGO_LW()).appendChild(sim_table)
        #xmldoc_mdc.appendChild(ligolw.LIGO_LW()).appendChild(sngl_table)
        ligolw_utils.write_filename(xmldoc_mdc, "mdc.xml", verbose = True)
        
        # build and write PSD XML document
        xmldoc_psd = lal.series.make_psd_xmldoc(PSDs,None)
        ligolw_utils.write_filename(xmldoc_psd, "psd.xml", verbose = True)
