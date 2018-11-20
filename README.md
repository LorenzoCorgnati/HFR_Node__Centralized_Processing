# HFR_Node_Centralized_Processing
Matlab scripts for the operational workflow of the European HFR Node. Tools for the centralized processing

These applications are written in Matlab language and they are based on HFR_Progs_2_1_2 and M_Map toolboxes, and the architecture of the workflow is based on a MySQL database containing information about data and metadata. The applications are designed for High Frequency Radar (HFR) data management according to the European HFR node processing workflow, thus generating radial and total velocity files in netCDF format according to the European standard data and metadata model for near real time HFR current data.

The database is composed by the following tables:
- account_tb: it contains the general information about HFR providers and the HFR networks they manage.
- network_tb: it contains the general information about the HFR network producing the radial and total files. These information will be used for the metadata content of the netCDF files.
- station_tb: it contains the general information about the radar sites belonging to each HFR network producing the radial and total files. These information will be used for the metadata content of the netCDF files.
- radial_input_tb: it contains information about the radial files to be converted and combined into total files.
- radial_HFRnetCDF_tb: it contains information about the converted radial files.
- total_input_tb: it contains information about the total files to be converted.
- total_HFRnetCDF_tb: it contains information about the converted total files.

The applications are intended to:
- load radial files information onto the database in table radial_input_tb;
- load total files information onto the database in table total_input_tb;
- convert Codar native .tuv files and WERA native .cur_asc files for total currents into the European standard data and metadata model for near real time HFR current data;
- convert Codar native .ruv files and WERA native .crad_ascii files for radial currents into the European standard data and metadata model for near real time HFR current data and combine them for generating total current files according to the European standard data and metadata model for near real time HFR current data.

General information for the tables network_tb and station_tb are loaded onto the database via a webform to be filled by the data providers. The webform is available at http://150.145.136.36/index.php

All generated radial and total netCDF files are quality controlled according the the QC tests defined as standard for the European HFR node and for the data distribution on CMEMS-INSTAC and SeaDataNet platforms.

The whole workflow is intended to run automatically to continuously convert and combine near real time HFR data produced by data providers. The wrapper CP_EU_HFR_Node_Processor.m sets the provider username and lauches the input and processing applications within an infinite loop.

The applications CP_inputRUV2DB.m and CP_inputCradAScii2DB.m load radial files information onto the database in table radial_input_tb.

The applications CP_inputTUV2DB.m and CP_inputCurAsc2DB.m load total files information onto the database in table total_input_tb.

The application CP_HFR_Combiner.m converts Codar native .ruv files and WERA native .crad_ascii files for radial currents into the European standard data and metadata model for near real time HFR current data and combines them for generating total current files according to the European standard data and metadata model for near real time HFR current data.

The application CP_Total_Conversion.m converts Codar native .tuv files and WERA native .cur_asc files for total currents into the European standard data and metadata model for near real time HFR current data.


The required toolboxes are:
- HFR_Progs-2.1.2 (https://github.com/rowg/hfrprogs); 
- M_Map (https://www.eoas.ubc.ca/~rich/map.html); 
- GSHHS (http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/); 
- Nctoolbox-1.1.3 (https://github.com/nctoolbox/nctoolbox); 
- mysql-connector-java-5.1.17 driver (https://mvnrepository.com/artifact/mysql/mysql-connector-java/5.1.17); 
- Rdir (http://www.mathworks.com/matlabcentral/fileexchange/19550-recursive-directory-listing);
- uniqueStrCell (https://www.mathworks.com/matlabcentral/fileexchange/50476-unique-for-cell-array-of-string).


Author: Lorenzo Corgnati

Date: November 20, 2018

E-mail: lorenzo.corgnati@sp.ismar.cnr.it
