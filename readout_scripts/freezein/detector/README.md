1. freezein_IC86_time_extraction.ipynb

is needed for extracting timestamps from the log file of each device run
is done by Leonhard already and saved in json file where logfile and raw device data is stored as well

2. freezein_IC86_extract.ipynb

loads timestamps from json file explained above and extracts the times,charges from all files/frames in the run time window
some is in /data/user/leidensc/uprade_installation_freezein feel free to load it from there or run it again and saved it somewhere else

3. freezein_IC86_readout.ipynb

loads the extracted data and the timestamps
bins the data and plots it

4. show_dev_day.ipynb & overview_plots.ipynb

is just freezein_IC86_readout.ipynb in a python script so i can load it in the notebook and loop over several days, devices etc...
