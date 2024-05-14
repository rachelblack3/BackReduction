import numpy as np
import matplotlib.pyplot as plt
import global_use as gl

# for cdf reading
from spacepy import pycdf

# other numerical tools
import os
from datetime import timedelta
from datetime import datetime


# Day range

start_date,end_date,no_days = gl.what_dates(default = True)
freq_bins = 65
time_stamps = 14400
freq_lists = np.zeros((no_days, time_stamps,freq_bins))

m = 0
for single_day in (start_date + timedelta(n) for n in range(no_days)):
    
    day_files = gl.DataFiles(single_day)

    # String version of dates for filename  
    date_string,year,month,day =gl.get_date_string(single_day)
    date_params = {"year": year,
                    "month": month,
                    "day": day,
                    "single_day": single_day}
    
    
    # Getting survey file and accessing survey frequencies, epoch and magnitude

    survey_file = pycdf.CDF(day_files.survey_data)
    survey_data = gl.AccessSurveyAttrs(survey_file)

    survey_freq = survey_data.frequency
    survey_epoch = survey_data.epoch_convert()
    Btotal = survey_data.Bmagnitude
    Etotal = survey_data.Emagnitude

    for j in range(freq_bins):
        max_time = np.shape(Etotal)[0]
        for n in range(time_stamps):
            freq_lists[m,n,j] = Etotal[n,j]



# Define empty list for total E field array 

Etotal = np.zeros(wna_file['EuEu'].shape)

freq_lists[m,n,j] = Etotal[n,j]
m=m+1
print('Done another',m)
