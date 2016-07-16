
Instructions to run the complete simulation
===========================================

1. Generate a simulated galaxy catalog

nohup octave -q generate_galaxies.m

Produces galaxies.txt and galaxies.mat

2. Generate a catalog of mergers at defined R_merger rate [23.5 Myr^-3 per MW galaxy]
   [note: uses the Glade catalog for mergers with distance < 200 Mpc]

nohup octave -q associate_merger_to_gal.m

Produces mergers.txt and mergers.mat

3. Generate XML tables for further processing with Bayestar/little_hope

ipython generate_xml_tables.py mergers.txt

Produces mdc.xml and psd.xml

Generates merger time, inclination, phase at coalescence, polarization angle

Select events with at least one single-detector SNR > threshold [= 4]

Power spectral densities are taken from the "Prospects" paper
http://arxiv.org/abs/1304.0670

aLIGO (mid-stage, low-sensitivity curve -- 80 Mpc)
adV (early-stage, low-sensitivity curve -- 20 Mpc)

4. Run little_hope

bayestar_sim_to_tmpltbank mdc.xml -o templates.xml
bayestar_littlehope --detector H1 --detector L1 --min-triggers 2 --snr-threshold 4.0 --reference-psd psd.xml --template-bank templates.xml --waveform "TaylorF2threePointFivePN" mdc.xml -o coinc.xml
bayestar_localize_coincs --waveform "TaylorF2threePointFivePN" --f-low 30 coinc.xml
bayestar_plot_allsky 0.toa_phoa_snr.fits.gz --contour 90 --radec long lat -o skymap.png
