# -*- coding: utf-8 -*-

from math import log10
from scipy.interpolate import interp1d

import json


with open("conf_K.json", "r+") as f:
    conf = json.load(f)

### INPUTS ### ------------------------------

data_f = conf["data_f"]   # Data from Hydrus, pay attebntion to the columns! (The data must be in the right place in csv file)
dots   = conf["dots"]     # Data from the retention function to be used for the interpolation for the center flux.

### OUTPUTS ### ------------------------------
output_new = conf["output_new"]       # This file is the most interesting, uses the proposed method for the calculations of K(h)


q_f = conf["q_f"]     # File with flux (q) at each time, provided from the weight of the scale. Always readed, but not used for calculations when: use_q_data = False
use_q_data = conf["use_q_data"]        # True uses q/2, False uses q based on tensiometer

ths = conf["ths"]   # Estimated thS, for the retention function

min_pF = conf["min_pF"]  # Value of pF were water content is thS

delta2 = conf["delta"]   # novo delta para calculo de z para o gradiente -- Método New

height = conf["height"]  # height of the column

# Distance between the center to the upper and lower tensiometer
dist = conf["dist"]


############################### ------------------------- # read retention points and transform into a list

z1 = height * 0.5 - dist
z2 = height * 0.5 + dist
zc = height * 0.5
zb = 0.5 * (zc + z1)

#import sys
#sys.exit()


# Generating list with retention function points
dotz = []
with open(dots, "r") as f: 
    for i in f.readlines():
        i = i.replace("\n","")
        ii = i.split(",")
        try:
            dotz.append([float(ii[0]),log10(float(ii[0])),float(ii[1])])
        except:
            dotz.append([ii[0],"pF",ii[1]])
            dotz.append([10**min_pF, min_pF, ths])
dotz.append([10**6.8, 6.8, 0.0])  # oven dry
dotz.pop(0)

# Flux data, if the flux is calculated based on the weight
if use_q_data:
    q_data = []
    with open(q_f, "r") as f: 
        for i in f.readlines():
            i = i.replace("\n","")
            ii = i.split(",")
            try:
                q_data.append([float(ii[0]),float(ii[1])])
            except:
                pass

# function to interpolate the retention function dots
def interp(hx, dotz, min_pF):
    """
    dotz = list with [ h, pF , theta] in crescent order in h, with h positive

    """
    if hx > 0.0:
        hx = log10(hx)
    else:
        hx = min_pF
    
    # Encontrando o intervalor
    for iv in range(len(dotz)):
        if iv == 0:
            continue
        
        if hx >= dotz[iv-1][1] and hx < dotz[iv][1]:
            index = iv
    
    h1, h2, th1, th2 = dotz[index-1][1], dotz[index][1], dotz[index-1][2], dotz[index][2]
    theta = (th2-th1)*(hx-h1)/(h2-h1) + th1
    
    return theta



############################### ------------------------- # read data file  ----------------------

data = []
with open(data_f, "r") as f:
    for i in f.readlines():
        i = i.replace("\n","")
        ii = i.split(",")
        try:
            #         time,      h2(top tensio) ,   h1 (botton tensiometer)
            # Condition to fullfil empty data
            if ii[1] == '':
                top = 9999999
            else:
                top = ii[1]
                
            if ii[2] == '':
                bot = 9999999
            else:
                bot = ii[2]
                
            iii = [float(ii[0]), abs(float(top)), abs(float(bot))]
            data.append(iii)
        except:
            pass

#  function to interpolate flux (when q/2)
if use_q_data:
    q_t = []
    q_q = []
    for i in range(len(q_data)):
        if i == 0:
            continue
        q_t.append(q_data[i][0])
        q_q.append(q_data[i][1]*height)  # Se dividir por 5 o arquivo
    
    q_func = interp1d(q_t, q_q, fill_value='extrapolate')


###############################  ---------------------------- Data list ---------------------


# Stablishing minimum and maximum values for the tensiometers reasings
min_h1 = min(data, key=lambda x: x[2])[2]
min_h2 = min(data, key=lambda x: x[1])[1]

max_h1 = max(data, key=lambda x: x[2])[2]
max_h2 = max(data, key=lambda x: x[1])[1]

min_h = max([min_h2, min_h1])
max_h = min([max_h2, max_h1])

# Generating list with values: h_mid , h_c-dt , h_c+dt
new_dat_geo = []
for i in range(len(data)):
    if (data[i][1] <= max_h) or (data[i][2] <= max_h):
        
        # searching for the intermediate value, "look-up table"
        val = data[i][2]
        t_bot = data[i][0]
        flag = 0
        for ii in range(len(data)):
            if ii == 0:
                continue
            if (val >= data[ii-1][1]) and (val <= data[ii][1]):
                index = ii
                flag = 1
                break
        
        if flag == 1:
            
            h21 = data[index-1][1]
            h22 = data[index][1]
            t1  = data[index-1][0]
            t2  = data[index][0]
            tx  = (val-h21)*(t2-t1)/(h22-h21) + t1
            
            ave_geo = ((t_bot**2 + tx**2)/2)**0.5
            
            Dt2 = abs(tx**2 - t_bot**2)
            
            dt_geo = Dt2 * delta2 / abs(z2 - z1)
            
            
            bg_geo = ave_geo**2 + dt_geo
            bl_geo = ave_geo**2 - dt_geo
            
            geo_gt = bg_geo**0.5
            geo_lt = bl_geo**0.5
            
            #               tempo,   h,  +dt  ,  -dt
            new_temp_geo = [ave_geo, val, geo_gt, geo_lt]
            
            new_dat_geo.append(new_temp_geo)

for i in range(len(new_dat_geo)):
    if new_dat_geo[i][1] >= max_h:
        t_max = new_dat_geo[i][0]
        break

# time,  pot,  +dt  ,  -dt
gt_lis = []
lt_lis = []
h_mid_lis = []
time = []
for i in range(len(new_dat_geo)):
    gt_lis.append(new_dat_geo[i][2])
    lt_lis.append(new_dat_geo[i][3])
    h_mid_lis.append(new_dat_geo[i][1])
    time.append(new_dat_geo[i][0])


# Interpolation functions of tension in the center of the sample
h_mid = interp1d(time, h_mid_lis,bounds_error=False,fill_value='extrapolate')
# Interpolation of the tension function "a little above" the center
f_gt = interp1d(gt_lis, h_mid_lis,bounds_error=False,fill_value='extrapolate')
# Interpolation of the tension function "a little below" the center
f_lt = interp1d(lt_lis, h_mid_lis,bounds_error=False,fill_value='extrapolate')
#  format:     x(time) , y(h)

# Geneerating Gradient list
grad = []
for i in range(len(new_dat_geo)):
    t = new_dat_geo[i][0]
    gg = (f_lt(t) - f_gt(t))/(2 * delta2)
    temp = [t,gg]
    grad.append(temp)
    

# Calculation of conductivity, at h_center 
data_out = []
verif = []
for i in range(len(data)-1):
    if data[i + 1][0] > t_max:
        break
    dt    = data[i + 1][0] - data[i][0]
    
    t1    = data[i][0]
    t2    = data[i + 1][0]
    
    h11   = data[i][2]     # tensiometer bottom (1), time i-1
    h12   = data[i + 1][2] # tensiometer bottom, time i
    th11  = interp(h11 , dotz , min_pF)  # water content of tensiometer bottom (1), time = i-1
    th12  = interp(h12 , dotz , min_pF)  # water content of tensiometer bottom (1), time = i
    
    h21   = data[i][1]      # tensiometer top (2), time i-1
    h22   = data[i + 1][1]  # tensiometer top, time i

    
    hmid1 = h_mid(t1)  # tension in the center of the sample, time = i-1
    hmid2 = h_mid(t2)  # tension in the center of the sample, time = i
    
    thmid1 = interp(hmid1 , dotz , min_pF)  # water content at the center of the soil sampel, time = i-1
    thmid2 = interp(hmid2 , dotz , min_pF)  # water content at the center of the soil sampel, time = i
    
    # Calculating displaced water ina  time-step
    l2 = (zc - zb) * thmid2 + zb * th12
    l1 = (zc - zb) * thmid1 + zb * th11
    d_th = abs(l2 - l1)
    
    q = d_th/dt
    
    t = data[i][0] + dt/2
    h_mid_temp = h_mid(t)
    
    gg = (f_lt(t) - f_gt(t))/(2 * delta2) - 1
    
    # Rewrite flux with flux from the balance, if desired
    if use_q_data:
        q = q_func(t)/2
    
    # Calculation of conductivity
    K = q/gg
    
    temp = [t, gg, q, K, h_mid_temp]
    data_out.append(temp.copy())
    


#########################   GENERATING FILE *.out  ----------------------------------------------------------------


with open(output_new, "w+") as f:
    f.write("time,grad,q,K,h_mid\n")
    for ii in range(len(data_out)):
        
        string = str(data_out[ii][0]) + ',' + str(data_out[ii][1]) + ',' + str(data_out[ii][2]) + ',' + str(data_out[ii][3]) + ',' + str(data_out[ii][4]) + '\n'
        f.write(string)








