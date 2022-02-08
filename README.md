# mosaic_mwr_data_desc
Visualization and creation of products derived from microwave radiometers during MOSAiC

Contact: a.walbroel@uni-koeln.de (MSc. Andreas Walbröl @University of Cologne, Institute for Geophysics and Meteorology)

Required packages (version): Python (3.8.10) 
tensorflow (2.5.0), keras (2.5.0), numpy (1.17.4 and 1.19.5 (latter for NN retrieval)), sklearn (0.24.2), netCDF4 (1.5.3 and 1.5.7 (latter for NN retrieval)), matplotlib (3.4.3), and xarray (0.18.2).

!!!!!
The scripts are linked to the publication (DOI_PUBLICATION). All scripts (except for data_tools.py, met_tools.py, import_data.py, my_classes.py) require to change paths. Paths can be found near the top of scripts (below the functions). 
!!!!!

Following data sets are linked to these scripts:
Retrieval training:
  https://doi.org/10.5281/zenodo.5846394
microwave radiometer data:
  DOI_MIRAC_P_TB
  DOI_MIRAC_P_L2
  DOI_HATPRO_L1
  DOI_HATPRO_L2
Polarstern track data:
  https://doi.org/10.1594/PANGAEA.924668
  https://doi.org/10.1594/PANGAEA.924674
  https://doi.org/10.1594/PANGAEA.924681
  https://doi.org/10.1594/PANGAEA.926829
  https://doi.org/10.1594/PANGAEA.926910
Polarstern radiosonde data:
  https://doi.org/10.1594/PANGAEA.928656

Listing the scripts and their purposes:
NN_retrieval_miracp.py: Training and application of a Neural Network retrieval to derive integrated water vapour from the microwave radiometer MiRAC-P. Training data consists of ERA-I reanalysis (https://doi.org/10.5281/zenodo.5846394). Furthermore, the brightness temperature data measured by MiRAC-P (DOI_MIRAC_P_TB) during the MOSAiC expedition (https://doi.org/10.25923/9g3v-xh92) are required. When running the code with "python3 NN_retrieval_miracp.py 'finalized'" the DOI_MIRAC_P_L2 data is created. When running it with "python3 NN_retrieval_miracp.py '20_runs'" metrics mentioned in the DOI_PUBLICATION can be found.

data_tools.py: Module containing data analysis routines called by other scripts.

import_data.py: Module containing various importer routines called by other scripts.

met_tools.py: Module containing meteorological computations (humidity conversion, ...) called by other scripts

my_classes.py: Classes called by other scripts

case_study_overview_mwr_radiosonde.py: Script to generate Figure 6 of DOI_PUBLICATION. Uses DOI_HATPRO_L2, DOI_MIRAC_P_L2, netCDF version of https://doi.org/10.1594/PANGAEA.928656 (converted with PANGAEA_tab_to_nc.py).

mwr_pro_output_add_geoinfo.py: Adding Polarstern track data (https://doi.org/10.1594/PANGAEA.924668, https://doi.org/10.1594/PANGAEA.924674, https://doi.org/10.1594/PANGAEA.924681, https://doi.org/10.1594/PANGAEA.926829, https://doi.org/10.1594/PANGAEA.926910 <<-- as netCDF versions, converted with PANGAEA_tab_to_nc.py) to HATPRO files (DOI_HATPRO_L1, DOI_HATPRO_L2).

PANGAEA_tab_to_nc.py: Script to convert PANGAEA radiosonde (https://doi.org/10.1594/PANGAEA.928656) and Polarstern track data (see above) to netCDF format (.nc).

plot_mwr_level_2a_radiosonde.py: Script to generate Figures 1, 2, and 3 of DOI_PUBLICATION. Uses DOI_HATPRO_L2, DOI_MIRAC_P_L2, https://doi.org/10.1594/PANGAEA.928656 (netCDF version from PANGAEA_tab_to_nc.py).

plot_mwr_level_2bc_radiosonde.py: Script to generate Figures 4, and 5 of DOI_PUBLICATION. Uses DOI_HATPRO_L2 and https://doi.org/10.1594/PANGAEA.928656 (netCDF version from PANGAEA_tab_to_nc.py).
