#%% package import
# package import
# ==========================================
import pyEW
import numpy as np
import pandas as pd
import matplotlib

import importlib as imp
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import warnings
import json
from SALib.sample import morris as morris_sample
from SALib.analyze import morris as morris_analyze
import datetime
import re
from scipy.optimize import differential_evolution
import validation as vd
import numpy as np

#%% Initial setting
# Initial setting
# ==========================================
#Experiment time
start_time = '2022-05-23 00:00:00'
end_time = '2022-09-19 23:59:59'

#units
conv_mol = 1e6 # Conversion from moles to micromols 
conv_Al = 1e3 # Conversion for Al species from mols to nanomols

#Soil Biogeochemistry
[MM_Mg, MM_Ca, MM_Na, MM_K, MM_Si, MM_carbon, MM_Anions, MM_Al,MM_Ni]=pyEW.MM(conv_mol) 
MM_conv = {'Ca': MM_Ca, 'Mg': MM_Mg, 'K': MM_K, 'Na': MM_Na} # molar masses [g/µmol]
ions = ['Ca', 'Mg', 'K', 'Na']
Group = ['Control - Biochar - Clay','Control - Clay','Olivine - Biochar - Clay','Olivine - Clay','Wollastonite - Biochar - Clay','Wollastonite - Clay']
Group_plot = [g.replace('Olivine', 'Dunite') for g in Group]
code_map = {
    'C': 'Control - Clay',
    'BC': 'Control - Biochar - Clay',
    'OC': 'Olivine - Clay',
    'OBC': 'Olivine - Biochar - Clay',
    'WC': 'Wollastonite - Clay',
    'WBC': 'Wollastonite - Biochar - Clay'
}

# ==========================================
#Initial setting
# ==========================================
#water balance
keyword_wb = 1 # 1 = varying soil moisture
#background inputs of cations and anions
keyword_add = 0 # 0 = no addition

#IBC setting
soil = "clay" 
pH_in = 6.76
rho_bulk = 1.54*1e6 #soil dry mass (g/m3)
Zr = 0.2 # [m]: soil depth
s_in = 0.45 #relative soil moisture, dont have the data, use the data from lab experiment
h_container = 0.85 # [m]
A_container = 1.16*0.96 # [m2]
V_container = h_container*A_container # [m3]
M_soil = 342.2 #kg soil

#initial pCO2
CO2_atm = pyEW.CO2_atm(conv_mol)
CO2_air_in = 50*pyEW.CO2_atm(conv_mol) # [mol-conv/l] 
ratio_aut_het = 0.5 # [-]: ratio of autotrophic to heterotrophic respiration 

#Effective CEC (measured)
CEC_tot = 28.2 # [mmol_c / 100 g dry-soil]
CEC_tot = CEC_tot*1e-5*rho_bulk*Zr*conv_mol # [mol_c]
#CEC fractions
f_Ca_in = 0.925
f_Mg_in = 0.124
f_K_in = 0.020
f_Na_in = 0.007
f_Al_in =  0.005
n_sand = 0.35  # sand porosity, from pyEW.soil_const('sand')
Si_in = 1093.83*1e-9*rho_bulk/MM_Si  # [micromol/l]  # [micromol/l]
Ni_react_init = 22e-6 / 59  # mg/kg → mol/kg
#base saturation
f_base_CEC= f_Ca_in + f_Mg_in + f_K_in + f_Na_in #c:/Users/11734/Desktop/Thesis/Data/field experiment/0.1 M BaCl2 extraction - 2022+2023.xlsx
f_acid_CEC= 1-f_base_CEC #f_Al_in + f_H_in
total_occupied = f_base_CEC + f_Al_in
if total_occupied > 1.0:
    scale_factor = 1.0 / total_occupied
    f_Ca_in *= scale_factor
    f_Mg_in *= scale_factor
    f_K_in  *= scale_factor
    f_Na_in *= scale_factor
    f_Al_in = f_Al_in*scale_factor
    f_H_in   = 0
else:
    f_H_in = max(0, 1.0 - f_base_CEC - f_Al_in)

#carbonates [micromol/m2]
CaCO3_in = 0.259 * (100.09/12.01) * 1.54e6 * 0.2 / 100.09 * 1e6 / 1e3  # 实测值
MgCO3_in = 0

#hydroclimatic
date_start = datetime.date(2022, 5, 23)
day1 = date_start.timetuple().tm_yday  # 143 #initial day for the simulation（based on the currunt month to calclulate the starting day)
latitude = 51.98*np.pi/180 # Renkum
albedo = 0.23 #use the data from lab experiment
altitude = 51 #[m] from Wikipedia
coastal = False 

#%% field Data 
# Field Data 
# ==========================================
#import real data
file_path_wind = r'C:/Users/11734/Desktop/Thesis/Data/wind speed from KNMI.txt' # real data from KNMI, station 275-Deelen
file_temp_and_rain = r'C:/Users/11734/Desktop/Thesis/Data/Weather station data/Sinderhoeve_Cleaned.xlsx' #Data from weather station Sinderhoeve
file_path_soil_temp =r"C:/Users/11734/Desktop/Thesis/Data/export-395e09b0-a120-4922-acc8-7ebb4551a844/Veenkampen_Soil.csv" #[°C]
file_path_vagetation = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/Microwave destruction - first maize harvest Sept 2022.xlsx"
file_path_watering = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/Overview water additions - Lysimeter experiment.xlsx"
file_path_PSD = r"c:/Users/11734/Desktop/Thesis/Data/field experiment/Particle Size Distribution with sieving above 2 mm and laser diffraction  - minerals and biochar.xlsx"
file_path_hourly = r"C:/Users/11734/Desktop/Thesis/Data/Hourly data.txt" #hourly temperature data in 2022 in Deelen
file_path_HNO3 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/0.43 M HNO3 extraction - 2022+2023.xlsx"
file_path_CaCl2 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/0.01 M CaCl2 extraction - 2022+2023.xlsx"
file_path_BaCl2 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/0.1 M BaCl2 extraction - 2022+2023.xlsx"
file_path_lysimeter_Ca = r'C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Ca.xlsx'
file_path_lysimeter_Mg = r'C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Mg.xlsx'
file_path_lysimeter_K = r'C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - K.xlsx'
file_path_lysimeter_Na = r'C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Na.xlsx'
file_path_Alk_pH = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/pH and alkalinity measurements rhizon samples - Lysimeter experiment.xlsx"
file_path_DIC = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/SFA data rhizons - Lysimeter experiment.xlsx"
file_path_SOC = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/LECO CN before and after LOI - 2022+2023.xlsx"
file_path_Si = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Si.xlsx"
file_path_Ni = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Ni.xlsx"

# ==========================================
#experimental data  
# ==========================================

Ni_soil_obs = {}
df_ni = pd.read_excel(file_path_CaCl2, sheet_name='Calc_figures')
df_ni.columns = df_ni.columns.astype(str).str.strip()
ni_col = next((c for c in df_ni.columns
               if 'ni' in c.lower() and 'kg' in c.lower()), None)
df_ni = df_ni[
    df_ni['Sampling date'].astype(str).str.contains('2022') &
    df_ni['Treatment'].astype(str).str.contains('Clay', case=False, na=False)
][['Treatment', 'Sampling date', ni_col]].copy()
df_ni['Treatment'] = df_ni['Treatment'].astype(str).str.strip()
df_ni[ni_col] = df_ni[ni_col] / 1000   # µg/kg → mg/kg
Ni_soil_obs['CaCl2'] = df_ni.groupby('Treatment')[ni_col].mean().to_dict()

#%% Real situation setting
# Real situation setting
# ==========================================
#windspeed
# ==========================================
header_line_index = 0
with open(file_path_wind, 'r') as f: #skip readme file 
    for i, line in enumerate(f):
        if line.strip().startswith('# STN,YYYYMMDD'):
            header_line_index = i
            break
df_wind = pd.read_csv(file_path_wind,delimiter=',',skiprows=header_line_index, skipinitialspace=True, parse_dates=['YYYYMMDD'])
df_wind.columns = df_wind.columns.str.replace('# ', '').str.strip()
df_wind.set_index('YYYYMMDD', inplace=True)
df_wind_subset = df_wind.loc[start_time:end_time]
t_end = len(df_wind_subset) # [d]: number of simulated days
dt = 1/(24*60) # [time resolution (d)]
t=np.arange(0,t_end,dt)
win_size = int(3/dt) #moving average time interval
steps_per_day = int(round(1/dt))
raw_wind = df_wind_subset['FG'].values/10 #FG = Daily mean wind speed in 0.1 m/s --changed into [m/s]
wind = np.repeat(raw_wind,steps_per_day)

# ==========================================
#temperature  
# ==========================================
header_line_index = 0
with open(file_path_hourly, 'r') as f: #skip readme file 
    for i, line in enumerate(f):
        if line.strip().startswith('# STN,YYYYMMDD'):
            header_line_index = i
            break
df_temp_KNMI = pd.read_csv(file_path_hourly, delimiter=',', skiprows=header_line_index, skipinitialspace=True)
df_temp_KNMI.columns = df_temp_KNMI.columns.str.replace('# ', '').str.strip()
dates = pd.to_datetime(df_temp_KNMI['YYYYMMDD'], format='%Y%m%d')
hours = pd.to_timedelta(df_temp_KNMI['HH'] - 1, unit='h')
df_temp_KNMI.index = dates + hours
df_hourly_2022 = df_temp_KNMI.loc['2022']

df_temp = pd.read_excel(file_temp_and_rain,sheet_name='2022') #lack data before May
df_temp.set_index('Date', inplace=True)
extended_start = pd.Timestamp(start_time) - pd.Timedelta(hours=1)
extended_end = pd.Timestamp(end_time) + pd.Timedelta(hours=1)
df_temp_subset = df_temp.loc[extended_start:extended_end]

temp_av = df_hourly_2022['T'].mean()/10 #2022 average temp

year_max = df_hourly_2022['T'].max()/10
year_min = df_hourly_2022['T'].min()/10
temp_ampl_yr = (year_max - year_min)/2 # [°C] yearly amplitude

daily_amplitudes = (df_hourly_2022['T'].resample('D').max()/10 - df_hourly_2022['T'].resample('D').min()/10)/2
temp_ampl_d = daily_amplitudes.mean()# [°C] daily amplitude

temp_series_1min = df_temp_subset['Temperature °C'].resample('1min').interpolate(method='linear')
temp_air = temp_series_1min.loc[start_time:end_time].values
 
df_1min = df_temp_subset.resample('1min').ffill()
temp_min = df_1min.groupby(df_1min.index.date)['Temperature °C'].transform('min').loc[start_time:end_time].values
temp_max = df_1min.groupby(df_1min.index.date)['Temperature °C'].transform('max').loc[start_time:end_time].values

df_raw_soil_temp = pd.read_csv(file_path_soil_temp,header=0, skiprows=[1])#TS_1_3_1 is the soil temperature under 20cm of grass 
df_raw_soil_temp['Timestamp'] = pd.to_datetime(df_raw_soil_temp['Timestamp'])
df_raw_soil_temp.set_index('Timestamp', inplace=True)
temp_soil = df_raw_soil_temp.loc[start_time:end_time,'TS_1_3_1']
temp_soil = temp_soil[~temp_soil.index.duplicated(keep='first')]
temp_soil = temp_soil.values

ET0 = pyEW.ET0(latitude,altitude,temp_air,temp_soil,temp_min,temp_max,wind,albedo,Zr,coastal,t_end,dt,day1)

# ==========================================
#vegetation
# ==========================================
T_v = t_end # [d] growth time
#In  addition, two seeds were sown in the middle of the plot to ensure eight plants per plot
df_drymass = pd.read_excel(file_path_vagetation,sheet_name='Massbl')
df_sand = df_drymass[df_drymass['Treatment'].str.contains('Clay', case=False, na=False)]
df_plants = df_sand.groupby(['Treatment', 'Plant nr'])['Dry weight biomass (g)'].sum().reset_index()
df_plants['Avg_Plant_Weight'] = df_plants.groupby('Treatment')['Dry weight biomass (g)'].transform('mean')
df_plants['k_v (g/m2)'] = df_plants['Avg_Plant_Weight']*8/A_container# [g/m2] carrying capacity (8 plants per lysimeter)
df_plants['v_in (g/m2)'] = 0.001*df_plants['k_v (g/m2)']  # [g/m2] seed biomass
df_plants = df_plants.groupby('Treatment', as_index=False).agg({'k_v (g/m2)': 'first',
                                                                'v_in (g/m2)': 'first'})
df_plants = df_plants[~df_plants['Treatment'].str.contains('k-feldspar', case=False, na=False)].reset_index(drop=True)
RAI = 0.75 #[m2/m2] root area index
root_d = 0.3*1e-3 #[m] average root diameter
t0_v = 0 # [d] starting day of growth

# ==========================================
#save parameters
# ==========================================
result = {}
for index, row in df_plants.iterrows():
    v = pyEW.veg(
        v_in = row['v_in (g/m2)'], 
        T_v  = T_v,
        k_v  = row['k_v (g/m2)'], 
        t0_v = t0_v,                 
        temp_soil = temp_soil,            
        dt   = dt
    )
    result[row['Treatment']] = {'v': v}


# ==========================================
#rain [m]
# ==========================================
rain = df_temp_subset['Rain mm'].resample('1min').interpolate(method='linear') 
rain = rain.loc[start_time:end_time].values/1000 # from mm to m, per minute and change into array

#add watering events
df_raw_watering = pd.read_excel(file_path_watering)
df_raw_watering['Date'] = pd.to_datetime(df_raw_watering['Date'], errors='coerce')
start_dt = pd.to_datetime(start_time).normalize()
end_dt   = pd.to_datetime(end_time).normalize()
mask = (df_raw_watering['Date'] >= start_time) & (df_raw_watering['Date'] <= end_time)
watering = df_raw_watering.loc[mask]
start_dt = pd.to_datetime(start_time).normalize()
end_dt   = pd.to_datetime(end_time).normalize()
dates = pd.date_range(start=start_dt, end=end_dt, freq='D')
watering = watering.reindex(dates, fill_value=0)
watering['Rain_event'] = watering['Water added (L)']*1e-3/A_container # m/day
watering = watering['Rain_event'].values
# build full daily index from simulation start to end
for i in range(0, len(rain)):
    day_index = int(i/(24*60))
    rain[i] = rain[i] + watering[day_index]/(24*60) # m/min

# ==========================================
#moisture balance (I, Q [m], L, T, E [m/d])
# ==========================================
for index, row in df_plants.iterrows():
    [s, s_w, s_i, I, L, T, E, Q, Irr, n] = pyEW.moisture_balance(
        rain = rain, 
        Zr   = Zr, 
        soil = soil, 
        ET0  = ET0, 
        v    = result[row['Treatment']]['v'], 
        k_v  = row['k_v (g/m2)'], 
        keyword_wb = keyword_wb,
        s_in = s_in,
        t_end= t_end,
        dt   = dt)
    result[row['Treatment']].update({'s': s, 's_w': s_w, 's_i': s_i, 'I': I, 'L': rain, 'T': T, 'E': E, 'Q': Q, 'Irr': Irr,'n': n})

# ==========================================
# Carbon system
# ==========================================

#Initial organic carbon
ADD = 0 # added dry litter [gOC/(m2*d)]
SOC_in = rho_bulk*13.9/1000 #[gOC/m3] initial sand OC content

#SOC balance and respiration
for index, row in df_plants.iterrows():
    treatment = row['Treatment']
    if 'Biochar' in treatment:
        Biochar_Add = rho_bulk*6.5/1000 # gOC/m3 soil added with biochar application
    else:
        Biochar_Add = 0

    [SOC_total, r_het, r_aut, D] = pyEW.respiration(
        ADD = ADD, 
        SOC_in = SOC_in, 
        CO2_air_in = CO2_air_in, 
        ratio_aut_het = ratio_aut_het, 
        soil = soil, 
        s = result[treatment]['s'], 
        v    = result[treatment]['v'], 
        k_v  = row['k_v (g/m2)'], 
        Zr = Zr, 
        temp_soil = temp_soil,
        dt = dt,
        conv_mol = conv_mol,
        Biochar_Add = Biochar_Add,
        k_dec_factor = 0.2)
    result[treatment].update({'SOC': SOC_total, 'r_het': r_het, 'r_aut': r_aut, 'D': D})

# ==========================================
# Particle size distribution
# ==========================================
data_exp_CDF = pd.read_excel(file_path_PSD,sheet_name='Calc',header=0)
data_exp_PSD = pd.read_excel(file_path_PSD,sheet_name='Sheet1',header=0)
size_row_index = 1
start_colomn_index = 21 #start from colomn20-'Result Between User Sizes (Sizes in um)'
num_bin = 101

Dunite_PSD = data_exp_PSD[data_exp_PSD['Sample Name']=='1_olivine Greensand'].values[0, start_colomn_index : start_colomn_index + num_bin].flatten()
Wollastonite_PSD = data_exp_PSD[data_exp_PSD['Sample Name']=='2_Wollastonite'].values[0, start_colomn_index : start_colomn_index + num_bin].flatten()
d_in = np.array(data_exp_PSD.iloc[0, start_colomn_index : start_colomn_index + num_bin].values.flatten())*1e-6 #[m]: particle diameter
psd_perc_in_wollastonite = np.array(Wollastonite_PSD/np.sum(Wollastonite_PSD)) #%
psd_perc_in_dunite = np.array(Dunite_PSD/np.sum(Dunite_PSD)) 
D50_dunite = 40.7e-6 #[m]: average diameter
D50_wollastonite = 108e-6 #[m]: average diameter

for i in Group:
    if 'Olivine' in i:
        result[i].update({'PSD_perc': Dunite_PSD, 'psd_perc_in': psd_perc_in_dunite, 'd50': D50_dunite,})
    elif 'Wollastonite' in i:
        result[i].update({'PSD_perc': Wollastonite_PSD, 'psd_perc_in': psd_perc_in_wollastonite, 'd50': D50_wollastonite})  
    elif 'Control' in i:
        result[i].update({'PSD_perc': np.array([0]),'psd_perc_in': np.array([0]),'d50': 0})

# ==========================================
#CEC
# ==========================================

for i in Group:
    result[i]['f_CEC_in'] = np.array([
    f_Ca_in,
    f_Mg_in,
    f_K_in,
    f_Na_in,
    f_Al_in,
    f_H_in
]) 
if abs(sum(result[i]['f_CEC_in'])-1) > 1e-3:
    raise ValueError("Sum of fractions must be 1")


for i in Group:
    [conc_in, K_CEC] = pyEW.f_CEC_to_conc(
        f_CEC_in = result[i]['f_CEC_in'], 
        pH_in = pH_in, 
        soil =  soil, 
        conv_mol = conv_mol,
        conv_Al = conv_Al)
    result[i].update({'conc_in': conc_in, 'K_CEC': K_CEC})

if abs(sum(result[i]['f_CEC_in'])-1) > 1e-6:
    raise ValueError("Sum of fractions must be 1")

def linear_CEC(CEC_start, CEC_end, n_steps):
    return np.linspace(CEC_start, CEC_end, n_steps)

CEC_endpoints = {
    'Control - Clay': {
        'start': CEC_tot,      
        'end': 29.59*1e-5*rho_bulk*Zr*conv_mol         
    },
    'Control - Biochar - Clay': {
        'start': CEC_tot,
        'end': 29.45*1e-5*rho_bulk*Zr*conv_mol           
    },
    'Olivine - Clay': {
        'start': CEC_tot,
        'end': 27.43*1e-5*rho_bulk*Zr*conv_mol    
    },
    'Olivine - Biochar - Clay': {
        'start': CEC_tot,
        'end': 29.73*1e-5*rho_bulk*Zr*conv_mol    
    },
    'Wollastonite - Clay': {
        'start': CEC_tot,
        'end': 30.56*1e-5*rho_bulk*Zr*conv_mol    
    },
    'Wollastonite - Biochar - Clay': {
        'start': CEC_tot,
        'end': 31.02*1e-5*rho_bulk*Zr*conv_mol    
    }
}

OM_pct = SOC_in / rho_bulk * 1000 / 1.724  # [%] SOC(g/kg) ÷ 1.724 = OM(%)

print('Initial soil setting finished.')
#%% Calibration
diss_f_dict = {
    'Olivine - Clay':                0.2216*0.1, #1.3%
    'Olivine - Biochar - Clay':      0.1759*0.1, #0
    'Wollastonite - Clay':           0.9353*0.1,#18.5%
    'Wollastonite - Biochar - Clay': 0.9353*0.1,#18.4%
    'Control - Clay':                0.0,
    'Control - Biochar - Clay':      0.0,
}

diss_f_n_dict = {
    'Olivine - Clay':                0.0,   # to be calibrated
    'Olivine - Biochar - Clay':      0.0,   # to be calibrated
    'Wollastonite - Clay':           0.0,   # to be calibrated
    'Wollastonite - Biochar - Clay': 0.0,   # to be calibrated
    'Control - Clay':                0.0,
    'Control - Biochar - Clay':      0.0,
}

k_Si_dict = {
    'Olivine - Clay':                0.03044, #0%
    'Olivine - Biochar - Clay':      0.03716, #0%
    'Wollastonite - Clay':           0.08201, #0%
    'Wollastonite - Biochar - Clay': 0.06086, #0.0%
    'Control - Clay':                0.07500, #71.7%,
    'Control - Biochar - Clay':      0.07500, #68.7%
}

print("Calibration complete. Using pre-calibrated values.")

#%% Run model
# ==========================================
def run_all_treatments():
    data = {}
    for index, row in df_plants.iterrows():
        i = row['Treatment']
        if i not in Group:
            continue

        if 'Olivine' in i:
            mineral   = ['forsterite','enstatite','fayalite','lizardite','clinochlore']
            rock_f_in = np.array([0.54, 0.02, 0.03, 0.24, 0.06])
            SSA_in    = 3.8
            M_rock_in = 26 * rho_bulk * Zr / 1000
            psd_perc  = psd_perc_in_dunite
            d_in_use  = d_in
        elif 'Wollastonite' in i:
            mineral   = ['wollastonite','diopside','albite','alkali_feldspar']
            rock_f_in = np.array([0.45, 0.13, 0.09, 0.15])
            SSA_in    = 0.9
            M_rock_in = 26 * rho_bulk * Zr / 1000
            psd_perc  = psd_perc_in_wollastonite
            d_in_use  = d_in
        else:
            mineral   = ['wollastonite']
            rock_f_in = np.array([0])
            SSA_in    = 0
            M_rock_in = 0
            psd_perc  = np.array([0])
            d_in_use  = np.array([0])

        data[i] = pyEW.biogeochem_balance(
            n           = result[i]['n'],
            s           = result[i]['s'],
            L           = result[i]['L'],
            T           = result[i]['T'],
            I           = result[i]['I'],
            v           = result[i]['v'],
            k_v         = row['k_v (g/m2)'],
            RAI         = RAI,
            root_d      = root_d,
            Zr          = Zr,
            r_het       = result[i]['r_het'],
            r_aut       = result[i]['r_aut'],
            D           = result[i]['D'],
            temp_soil   = temp_soil,
            pH_in       = pH_in,
            conc_in     = result[i]['conc_in'],
            f_CEC_in    = result[i]['f_CEC_in'],
            K_CEC       = result[i]['K_CEC'],
            CEC_tot     = linear_CEC(CEC_endpoints[i]['start'],
                                     CEC_endpoints[i]['end'], len(t)),
            Si_in       = Si_in,
            CaCO3_in    = CaCO3_in,
            MgCO3_in    = MgCO3_in,
            M_rock_in   = M_rock_in,
            t_app       = 0,
            mineral     = mineral,
            rock_f_in   = rock_f_in,
            d_in        = d_in_use,
            psd_perc_in = psd_perc,
            SSA_in      = SSA_in,
            diss_f      = diss_f_dict[i],
            dt          = dt,
            conv_Al     = conv_Al,
            conv_mol    = conv_mol,
            keyword_add = keyword_add,
            k_prec_Si   = k_Si_dict[i],
            keyword_Si_precip = 1,
            rho_bulk       = rho_bulk,
            Ni_react_init  = Ni_react_init,
            OM_pct         = OM_pct,
        )
        print(f'  {i} finished.')
    return data

print('Running model...')
biogeochem_out = run_all_treatments()
print('Model finished.\n')

#%% Unit conversion - Ni
# ==========================================
for i in Group:
    data = biogeochem_out[i]
    result[i]['Ni_soil'] = pyEW.mov_avg(
        data['Ni_tot'] * MM_Ni * 1000 / (M_soil / A_container), win_size)  # [mg/kg]
    result[i]['Ni_pw']   = pyEW.mov_avg(data['Ni'], win_size)               # [µmol/L]

print('Unit conversion finished.')

#%% Ni soil figure
# ==========================================
# 读取 Ni 土壤观测数据（HNO3 & CaCl2 提取）
# ==========================================
Ni_soil_obs = {}
for ext_name, fpath in [('HNO3', file_path_HNO3), ('CaCl2', file_path_CaCl2)]:
    df = pd.read_excel(fpath, sheet_name='Calc_figures')
    df.columns = df.columns.astype(str).str.strip()
    ni_col = next((c for c in df.columns
                   if 'ni' in c.lower() and 'kg' in c.lower()), None)
    if ni_col is None:
        print(f'Warning: {ext_name} 中没有找到 Ni 列，跳过')
        continue
    df = df[
        df['Sampling date'].astype(str).str.contains('2022') &
        df['Treatment'].astype(str).str.contains('Clay', case=False, na=False)
    ][['Treatment', 'Sampling date', ni_col]].copy()
    df['Treatment'] = df['Treatment'].astype(str).str.strip()
    if 'µg' in ni_col or 'ug' in ni_col.lower():
        df[ni_col] = df[ni_col] / 1000
    Ni_soil_obs[ext_name] = df.groupby('Treatment')[ni_col].mean().to_dict()

# ==========================================
# 作图
# ==========================================
colors = ['#8FA5C7','#4C72B0','#EDB68A','#DD8452','#F77DE8','#F00BD9']
color  = {g: c for g, c in zip(Group, colors)}
markers_method = {'HNO3': 'o', 'CaCl2': '^'}

fig, axes = plt.subplots(1, 3, figsize=(21, 7))
ax_con, ax_woll, ax_ol = axes

for i in Group:
    if 'Control' in i:
        ax = ax_con
    elif 'Wollastonite' in i:
        ax = ax_woll
    elif 'Olivine' in i:
        ax = ax_ol
    else:
        continue

    ls = '-' if 'Biochar' in i else '--'
    c  = color[i]

    # 模型线
    ax.plot(t, result[i]['Ni_soil'], color=c, linewidth=2,
            linestyle=ls, alpha=0.8, label=i)

    # 观测点（实验开始 day=0，实验结束 day=t[-1]）
    for ext_name, marker in markers_method.items():
        if ext_name not in Ni_soil_obs:
            continue
        fc = c if 'Biochar' in i else 'none'
        ec = 'white' if 'Biochar' in i else c
        lw_m = 1 if 'Biochar' in i else 1.5

        if 'Start soil sand' in Ni_soil_obs[ext_name]:
            ax.scatter(0, Ni_soil_obs[ext_name]['Start soil sand'],
                       color=c, marker=marker, s=80,
                       facecolors='none', edgecolors=c,
                       linewidths=1.5, zorder=5)
        if i in Ni_soil_obs[ext_name]:
            ax.scatter(t[-1], Ni_soil_obs[ext_name][i],
                       color=c, marker=marker, s=80,
                       facecolors=fc, edgecolors=ec,
                       linewidths=lw_m, zorder=5)

# 坐标轴设置
titles = ["Control (Clay)", "Wollastonite (Clay)", "Olivine/Dunite (Clay)"]
for ax, title in zip(axes, titles):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_ylim(bottom=0)
axes[0].set_ylabel("Ni in soil (mg/kg)", fontsize=12)

# 图例
legend_model = [
    mlines.Line2D([], [], color=color[g], linewidth=2,
                  linestyle='-' if 'Biochar' in g else '--', label=g)
    for g in Group
]
legend_obs = [
    mlines.Line2D([], [], color='black', marker=markers_method[m],
                  linestyle='None', markersize=8, label=f'{m} observed')
    for m in markers_method
]

axes[1].legend(handles=legend_model + legend_obs + legend_biochar,
               loc='upper center', bbox_to_anchor=(0.5, -0.15),
               ncol=4, frameon=False, fontsize=10)
plt.tight_layout()
plt.savefig('Ni_soil.png', dpi=150, bbox_inches='tight')
plt.show()
# %%
import re

file_path = r'C:/Users/11734/Desktop/Thesis/SMEW-Thesis - Ni/Clay - Thesis.py'
with open(file_path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

for idx, line in enumerate(lines, 1):
    if 'Si_in' in line and '=' in line:
        print(f"Line {idx}: {line.rstrip()}")
# %%
