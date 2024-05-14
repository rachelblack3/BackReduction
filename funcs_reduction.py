import numpy as np 
import matplotlib.pyplot as plt
from datetime import datetime,date,timedelta
import os
#os.environ["CDF_LIB"] = "/users/rablack75/"
from spacepy import pycdf
import glob
from scipy.optimize import curve_fit
import global_use as gu

def gaussian_fit(x, a, b, c, offset):
    return a * np.exp(-0.5 * np.power((x-b) / c, 2.0)) + offset



# want to write a code that searches through all WFR electric field wave power files
# between the dates of 30th sept 2012 - 16th feb 2015

# finds the peak power at each frequency bin 
# plot peak power against frequency 

# First, need generic path name to VA probe data 
# with abillity to cycle through the dates 


va_data = '/data/spacecast/satellite/RBSP/emfisis/data/RBSP-A/L2'
file_name = 'rbsp-a_WFR-spectral-matrix-diagonal_emfisis-L2_'


#Initializing the first date and time

start_date = date(year=2012, month=9, day=30)

#Initializing the end date and time

end_date = date(year=2014, month=9, day=16)

#Calculating no of days in range 

no_days = end_date - start_date
no_days = int(no_days.days)

start_date_file = 20120930
end_date_file = 20150216


freq_bins = 65
time_stamps = 14400
freq_lists = np.zeros((10, time_stamps,freq_bins))


# Make data frame with 65 rows, and a value for each day
m=0
for single_date in (start_date + timedelta(n) for n in range(10)):
    
    date_string= str(single_date.strftime("%Y%m%d"))

    if (single_date.day <10):
        day = "0"+str(single_date.day)
    else:
        day = str(single_date.day)


    if (single_date.month<10):
        month = "0"+str(single_date.month)
    else:
        month = str(single_date.month)
    

    year = str(single_date.year)
    
    path = os.path.join(va_data, year, month, day,file_name + date_string + "_v*.cdf")
    # find the latest version
    path = glob.glob(path)[-1]
    

    wna_file= pycdf.CDF(path)

    epoch_t = wna_file['Epoch']
   # for i in range(len(epoch_t)):
        #epoch_t[i]=dt.strftime(epoch_t[i],'%Y-%m-%d %H-%M-%S.%f')
       # epoch_t[i]=dt.strptime(epoch_t[i],'%Y-%m-%d %H-%M-%S.%f')
    
# Extract all electric field components

    Eu2 = wna_file['EuEu']
    Ev2 = wna_file['EvEv']
    Ew2 = wna_file['EwEw']

# Define empty list for total E field array 

    Etotal = np.zeros(wna_file['EuEu'].shape)

    freq_lists[m,n,j] = Etotal[n,j]
    m=m+1
    print('Done another',m)

# Now we need to bin each frequency list based upon power 
# (power vs counts)
freq_i=[]
freq_log=[]
freq_log_i=[]
sd_i=[]
N =[]
bw=[]
no_bin=[]
max_power=[]
sd_2=[]
sd_2_i=[]
pre_std=[]
fit_y=[]
x_i=[]
reg_y=[]
x_mod=[]


for i in range(65):
    print(i)
    freq_i.append(list(freq_lists[:,:,i].flatten()))
    freq_i[i] = list(filter(lambda x: x != -1., freq_i[i]))
    freq_i[i] = list(filter(lambda x: x != 0., freq_i[i]))
    

    for g in range(len(freq_i[i])):
        freq_log.append(np.log10(freq_i[i][g]))
    
    freq_log_i.append(freq_log)
    # log the data, then use Scott (1979) rule for bin width
    sd_i.append(np.std(freq_log_i[i]))
    N.append(len(freq_log_i[i]))
    bw.append(3.49*sd_i[i]/(N[i]**(1./3.)))
    no_bin.append(int((max(freq_log_i[i])-min(freq_log_i[i]))/bw[i]))

    # create histogram

    y= list(np.histogram(freq_log_i[i],bins=no_bin[i])[0])
    x= list(np.histogram(freq_log_i[i],bins=no_bin[i])[1])
    

    # find middle value of each bin 

    x_mid=[]
    for l in range(len(x)-1):
        
        x_mid.append(x[l]+(x[l+1]-x[l])/2)
    
    txt_file = "freq_no"+str(i)+".dat"
    data_save = np.column_stack((x_mid,y))
    np.savetxt(txt_file, data_save)

    max_y = int(y.index(max(y)))

    gauss_y = []
    gauss_x=[]

    for x_val in range(0,max_y+1,1):

        gauss_x.append(x_mid[x_val])
        gauss_y.append(y[x_val])
    
    for x_val in range(0,max_y,1):

        gauss_x.append(gauss_x[max_y]-gauss_x[max_y-x_val]+gauss_x[max_y])
        gauss_y.append(gauss_y[max_y-x_val])

    

    # initial parameter estimates from the data for gaussian fit
    a = max(gauss_x)
    b = x_mid[max_y]
    c = 1.0 # my guess from the equation
    offset = min(gauss_y)
    initialParameters = np.array([a, b, c, offset])

    # curve fit the test data
    fittedParameters, pcov = curve_fit(gf, gauss_x, gauss_y, initialParameters)

    yfit = gf(gauss_x, *fittedParameters) 
    #absError = yfit - y
    sd_2.append(fittedParameters[1]+abs(2*fittedParameters[2]))
   
    fit_y.append(yfit)
    reg_y.append(y)
    x_mod.append(gauss_x)
    x_i.append(x_mid)
    
    max_power.append(fittedParameters[1])
    
    plt.figure()
    plt.plot(x_mod[i], fit_y[i], '-')
    plt.plot(x_i[i], reg_y[i], '-')
    plt.axvline(sd_2[i],color='black',linestyle='dashed')
    plt.xlabel(r'$log10(V^2/m^2/Hz)$')
    plt.ylabel('Counts')
    title="Frequency band"+str(i)
    name = title+".png"
    plt.title(title)
    #plt.gca().set_xscale("log")
    plt.savefig(name)
    

data = np.column_stack((sd_2,max_power))
np.savetxt('hockey_plot.dat', data)
np.savetxt('freq.dat',wna_file['WFR_frequencies'][0])
#plt.hist(freq_log_i[0], bins=no_bin[0])
plt.plot(x_mod[32], fit_y[32], '-')
plt.plot(x_i[32], reg_y[32], '-')
plt.axvline(sd_2[32],color='black',linestyle='dashed')
plt.xlabel(r'$log10(V^2/m^2/Hz)$')
plt.ylabel('Counts')
plt.title('Frequency middle band - 2.8 x 10^2 Hz ')
#plt.gca().set_xscale("log")
plt.savefig('histogram.png')
#plt.hist(freq_log_i[60], bins=no_bin[60])

#plt.savefig('histogram65.png')
#plt.plot(wna_file['WFR_frequencies'][0],max_power)
#plt.plot(wna_file['WFR_frequencies'][0],sd_2)
#plt.gca().set_xscale("log")
#plt.savefig('hockey_stick.png')