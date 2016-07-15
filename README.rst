
Instructions to run the complete simulation
===========================================

1. Generate a simulated galaxy catalog

nohup octave -q generate_galaxies.m

Produces galaxies.txt

2. Generate a catalog of mergers at defined R_merger rate
   [note: uses the Glade catalog as well]

nohup octave -q associate_merger_to_gal.m

Produces mergers.txt

3. Generate XML tables for further processing with Bayestar/little_hope

ipython generate_xml_tables.py mergers.txt

Produces mdc.xml and psd.xml

4. Run little_hope

bayestar_sim_to_tmpltbank mdc.xml -o templates.xml
bayestar_littlehope --detector H1 --detector L1 --min-triggers 2 --snr-threshold 4.0 --reference-psd psds_2016-17.xml --template-bank templates.xml --waveform "TaylorF2threePointFivePN" mdc.xml -o coinc.xml
bayestar_localize_coincs --waveform "TaylorF2threePointFivePN" --f-low 30 coinc.xml
bayestar_plot_allsky 0.toa_phoa_snr.fits.gz --contour 90 --radec long lat -o skymap.png
