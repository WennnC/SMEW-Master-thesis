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
import gc

#%% Initial setting
# Initial setting
# ==========================================
#Experiment time
start_time = '2022-05-23 00:00:00'
end_time = '2024-05-06 23:59:59'

#units
conv_mol = 1e6 # Conversion from moles to micromols 
conv_Al = 1e3 # Conversion for Al species from mols to nanomols

#Soil Biogeochemistry
[MM_Mg, MM_Ca, MM_Na, MM_K, MM_Si, MM_carbon, MM_Anions, MM_Al]=pyEW.MM(conv_mol) 
MM_conv = {'Ca': MM_Ca, 'Mg': MM_Mg, 'K': MM_K, 'Na': MM_Na} # molar masses [g/µmol]
ions = ['Ca', 'Mg', 'K', 'Na']
Group = ['Control - Clay','Control - Biochar - Clay','Olivine - Clay','Olivine - Biochar - Clay','Wollastonite - Clay','Wollastonite - Biochar - Clay']
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
s_in = 0.5 #relative soil moisture, dont have the data, use the data from lab experiment
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
f_Ca_in = 0.859
f_Mg_in = 0.115
f_K_in = 0.019
f_Na_in = 0.007
f_Al_in =  0.001
Si_in = 2.53 / 28 * 1000 #[micromol/l]
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
    f_Al_in = f_Al_in*scale_factor/2
    f_H_in   = f_Al_in
else:
    f_H_in = max(0, 1.0 - f_base_CEC - f_Al_in)

#carbonates [micromol/m2]
CaCO3_in = 4.75e6
MgCO3_in = 0

#hydroclimatic
date_start = datetime.date(2022, 5, 23)
day1 = date_start.timetuple().tm_yday  # 143 #initial day for the simulation（based on the currunt month to calclulate the starting day)
latitude = 51.98*np.pi/180 # Renkum
albedo = 0.23 #use the data from lab experiment
altitude = 51 #[m] from Wikipedia
coastal = False 

#%% Import file path
# Import file path 
# ==========================================
#import real data
file_path_wind = r'C:/Users/11734/Desktop/Thesis/Data/wind speed from KNMI.txt' # real data from KNMI, station 275-Deelen
file_temp_and_rain = r'C:/Users/11734/Desktop/Thesis/Data/Weather station data/Sinderhoeve_Cleaned.xlsx' #Data from weather station Sinderhoeve
file_path_soil_temp =r"C:/Users/11734/Desktop/Thesis/Data/export-395e09b0-a120-4922-acc8-7ebb4551a844/Veenkampen_Soil.csv" #[°C]
file_path_vegetation1 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/Microwave destruction - first maize harvest Sept 2022.xlsx"
file_path_vegetation2 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/Microwave destruction - first winter rye harvest May 2023.xlsx"
file_path_vegetation3 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/Microwave destruction - second maize harvest Sept 2023.xlsx"
file_path_vegetation4 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/Dry weight second rye harvest - Lysimeter experiment - May 2024 - Emily te Pas.xlsx"
file_path_watering = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/Overview water additions - Lysimeter experiment.xlsx"
file_path_PSD = r"c:/Users/11734/Desktop/Thesis/Data/field experiment/Particle Size Distribution with sieving above 2 mm and laser diffraction  - minerals and biochar.xlsx"
file_path_hourly = r"C:/Users/11734/Desktop/Thesis/Data/Hourly data.txt" #hourly temperature data in Deelen
file_path_HNO3 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/0.43 M HNO3 extraction - 2022+2023.xlsx"
file_path_CaCl2 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/0.01 M CaCl2 extraction - 2022+2023.xlsx"
file_path_BaCl2 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/0.1 M BaCl2 extraction - 2022+2023.xlsx"
file_path_2024 = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/All lab data final soil samples PhD - Lysimeter experiment - May 2024.xlsx"
file_path_lysimeter_Ca = r'C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Ca.xlsx'
file_path_lysimeter_Mg = r'C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Mg.xlsx'
file_path_lysimeter_K = r'C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - K.xlsx'
file_path_lysimeter_Na = r'C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Na.xlsx'
file_path_Alk_pH = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/pH and alkalinity measurements rhizon samples - Lysimeter experiment.xlsx"
file_path_DIC = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/SFA data rhizons - Lysimeter experiment.xlsx"
file_path_SOC = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/LECO CN before and after LOI - 2022+2023.xlsx"
file_path_lysimeter_Si = r"C:/Users/11734/Desktop/Thesis/Data/field experiment/ICP data rhizons - Lysimeter experiment - Si.xlsx"

#%%# Field data  
# Import experimental data  
# ==========================================
#cation in soil with different measurement
#2022 & 2023
path_dict = {
    'HNO3': file_path_HNO3,   
    'CaCl2': file_path_CaCl2,
    'BaCl2': file_path_BaCl2
}
cols_hno3_bacl2 = ['Treatment', 'Sampling date', 'Ca (mg/kg)', 'Mg (mg/kg)', 'K (mg/kg)', 'Na (mg/kg)']
cols_cacl2 = ['Treatment', 'Sampling date', 'Mg (mg/kg)', 'K (mg/kg)', 'Na (mg/kg)'] #no Ca data in CaCl2 extraction
real_data_dict = {}
for name, path in path_dict.items():
    df = pd.read_excel(path, sheet_name='Calc_figures')
    df = df[df['Sampling date'].astype(str).str.contains('May') & df['Treatment'].astype(str).str.contains('Clay', case=False, na=False)]
    if name == 'CaCl2':
        df = df[cols_cacl2]
    else:
        df = df[cols_hno3_bacl2]
    df = df.groupby(['Treatment', 'Sampling date'], as_index=False).mean()
    df = df.set_index('Treatment')
    real_data_dict[name] = df

#2024
for sheet in path_dict.keys():
    df = pd.read_excel(file_path_2024, sheet_name=sheet)
    if sheet == 'BaCl2':
        # cmol+/kg → mg/kg
        for col in ['Ca (cmol+/kg)', 'Mg (cmol+/kg)', 'K (cmol+/kg)', 'Na (cmol+/kg)']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df['Ca (mg/kg)'] = df['Ca (cmol+/kg)'] / 2 * 40.08 * 10
        df['Mg (mg/kg)'] = df['Mg (cmol+/kg)'] / 2 * 24.31 * 10
        df['K (mg/kg)']  = df['K (cmol+/kg)']  / 1 * 39.10 * 10
        df['Na (mg/kg)'] = df['Na (cmol+/kg)'] / 1 * 22.99 * 10
    df = df[df['Treatment'].astype(str).str.contains('0-20') & df['Treatment.1'].astype(str).str.contains('Clay', case=False, na=False)]
    df = df.drop(columns=['Treatment'])
    df = df.rename(columns={'Treatment.1': 'Treatment'})
    df['Sampling date'] = 'May_2024'
    if sheet == 'CaCl2':
        df = df[cols_cacl2]
    else:
        df = df[cols_hno3_bacl2]
    df = df.set_index('Treatment')
    real_data_dict[sheet] = pd.concat([real_data_dict[sheet], df])

# ==========================================
#cation in porewater and leachate
start_date = pd.Timestamp(start_time)
end_date = pd.Timestamp(end_time)
processed_data = {}
for ion in ions: 
    file_path = globals()[f'file_path_lysimeter_{ion}']
    try:           
        df_meta = pd.read_excel(file_path, sheet_name=f'{ion}_mgL', header=None)
        days_row = df_meta.iloc[0, 1:] 
        dates_row = df_meta.iloc[2, 1:]
        valid_days = []
        for day, date_val in zip(days_row, dates_row):
            try:
                d_num = float(day)
                date_ts = pd.to_datetime(date_val, errors='coerce')
                if pd.notna(date_ts) and (start_date <= date_ts <= end_date):
                    valid_days.append(d_num)
            except Exception as e:
                    continue
    except Exception as e:
            print(f"  failed to read the date: {e}")
            continue
        
    # select time
    df_data_25 = pd.read_excel(file_path,sheet_name=f'{ion}_mmolL_25',header=1)
    df_data_12 = pd.read_excel(file_path,sheet_name=f'{ion}_mmolL_12',header=1)
    df_data_25 = df_data_25.set_index(df_data_25.columns[0])
    df_data_12 = df_data_12.set_index(df_data_12.columns[0])
    cols_to_keep = []
    for c in df_data_25.columns:
        try:
            if float(c) in valid_days: 
                cols_to_keep.append(c)
        except ValueError:
            continue
    df_data_25 = df_data_25[cols_to_keep]
    df_data_12 = df_data_12[cols_to_keep]


    #select row
    valid_codes = set(code_map.keys())
    rows_to_keep_12 = []
    rows_to_keep_25 = []
    for sample_name_25 in df_data_25.index:
        code = sample_name_25.strip().split('-')
        if any(i in valid_codes for i in code):
            rows_to_keep_25.append(sample_name_25)
    for sample_name_12 in df_data_12.index:
        code = sample_name_12.strip().split('-')
        if any(i in valid_codes for i in code):
            rows_to_keep_12.append(sample_name_12)
    df_data_25 = df_data_25.loc[rows_to_keep_25]
    df_data_12 = df_data_12.loc[rows_to_keep_12]
    processed_data[f'{ion}_25cm'] = df_data_25
    processed_data[f'{ion}_12cm'] = df_data_12

# ==========================================
#Alk
df_Alk_25 = pd.read_excel(file_path_Alk_pH,sheet_name='Alk_25',header=1)
df_Alk_12 = pd.read_excel(file_path_Alk_pH,sheet_name='Alk_12',header=1)
df_Alk_25 = df_Alk_25.set_index(df_Alk_25.columns[0])
df_Alk_12 = df_Alk_12.set_index(df_Alk_12.columns[0])
final_cols = [c for c in cols_to_keep if c in df_Alk_25.columns]
df_Alk_25 = df_Alk_25[final_cols]
final_cols = [c for c in cols_to_keep if c in df_Alk_12.columns]
df_Alk_12 = df_Alk_12[final_cols]

#select row
final_rows_25 = [r for r in rows_to_keep_25 if r in df_Alk_25.index]
df_Alk_25 = df_Alk_25.loc[final_rows_25]
final_rows_12 = [r for r in rows_to_keep_12 if r in df_Alk_12.index]
df_Alk_12 = df_Alk_12.loc[final_rows_12]
processed_data['Alk_25cm'] = df_Alk_25
processed_data['Alk_12cm'] = df_Alk_12

# ==========================================
#DIC
df_DIC = pd.read_excel(file_path_DIC,sheet_name='IC_mmolL',header=0)
df_DIC = df_DIC.set_index(df_DIC.columns[0])
final_cols = [c for c in cols_to_keep if c in df_DIC.columns]
df_DIC = df_DIC[final_cols]

#select row
rows_to_keep_25 = []
rows_to_keep_12 = []
for sample_name in df_DIC.index:
    if isinstance(sample_name, str):
        code = sample_name.strip().split('-')
        if any(i in valid_codes for i in code) and '25' in code:
            rows_to_keep_25.append(sample_name)
        elif any(i in valid_codes for i in code) and '12.5' in code:
            rows_to_keep_12.append(sample_name)
df_DIC_25 = df_DIC.loc[rows_to_keep_25]
processed_data['DIC_25cm'] = df_DIC_25
df_DIC_12 = df_DIC.loc[rows_to_keep_12]
processed_data['DIC_12cm'] = df_DIC_12

# ==========================================
#pH
df_pH = pd.read_excel(file_path_Alk_pH,sheet_name='pH_overview',header=0)
df_pH = df_pH.set_index(df_pH.columns[1])
final_cols = [c for c in cols_to_keep if c in df_pH.columns]
df_pH = df_pH[final_cols]

#select row
df_pH_25 = df_pH.iloc[19:25]
processed_data['pH_25cm'] = df_pH_25
df_pH_12 = df_pH.iloc[0:6]
processed_data['pH_12cm'] = df_pH_12

# ==========================================
#CEC 2022 & 2023
df_CEC = pd.read_excel(file_path_BaCl2, sheet_name='Calc_bs')
df_CEC = df_CEC.set_index(df_CEC.columns[2])
df_extracted = df_CEC[['Sampling date', 'BS% (4)', 'CEC (cmol+/kg)', 'Ca %', 'Mg %', 'K %', 'Na %']].copy()
df_extracted['CEC(mmol+/kg)'] = df_extracted['CEC (cmol+/kg)'] * 10
df_extracted['alk_CEC'] = (df_extracted['BS% (4)'] / 100) * df_extracted['CEC(mmol+/kg)'] * M_soil
df_extracted = df_extracted[df_extracted['Sampling date'].astype(str).str.contains('May')& df_extracted.index.astype(str).str.contains('Clay', case=False, na=False)]
df_extracted = df_extracted.groupby(['Treatment', 'Sampling date']).mean(numeric_only=True).reset_index()
df_extracted = df_extracted.set_index('Treatment')

# 2024
df_CEC_2024 = pd.read_excel(file_path_2024, sheet_name='BaCl2')
df_CEC_2024 = df_CEC_2024[df_CEC_2024['Treatment'].astype(str).str.contains('0-20') & df_CEC_2024['Treatment.1'].astype(str).str.contains('Clay', case=False, na=False)]
df_extracted_2024 = df_CEC_2024[['Treatment.1', 'CEC (cmol+kg)', 'BS %', 'Ca (cmol+/kg)', 'Mg (cmol+/kg)', 'K (cmol+/kg)', 'Na (cmol+/kg)']].copy()
df_extracted_2024 = df_extracted_2024.rename(columns={'Treatment.1': 'Treatment'})
df_extracted_2024['Sampling date'] = 'May_2024'
for ion in ions:
    df_extracted_2024[f'{ion} %'] = df_extracted_2024[f'{ion} (cmol+/kg)'] / df_extracted_2024['CEC (cmol+kg)'] * 100
df_extracted_2024['alk_CEC'] = (df_extracted_2024['BS %'] / 100) * df_extracted_2024['CEC (cmol+kg)'] * M_soil
df_extracted_2024['CEC(mmol+/kg)'] = df_extracted_2024['CEC (cmol+kg)'] * 10
df_extracted_2024 = df_extracted_2024.set_index('Treatment')
common_cols = df_extracted.columns.intersection(df_extracted_2024.columns).tolist()
df_extracted = pd.concat([df_extracted, df_extracted_2024])[common_cols]

# ==========================================
#SOC 2022 & 2023
df_SOC = pd.read_excel(file_path_SOC, sheet_name='Calc')
df_SOC = df_SOC[
    df_SOC['Sampling date'].astype(str).str.contains('May') &
    df_SOC['Treatment'].astype(str).str.contains('Clay', case=False, na=False)
]
df_SOC = df_SOC[['Treatment', 'Sampling date', 'SOC (g/kg)']].copy()
df_SOC = df_SOC.groupby(['Treatment', 'Sampling date']).mean(numeric_only=True).reset_index()
df_SOC = df_SOC.set_index('Treatment')

# 2024
df_SOC_2024 = pd.read_excel(file_path_2024, sheet_name='LECO')
df_SOC_2024 = df_SOC_2024[df_SOC_2024['Treatment'].astype(str).str.contains('0-20') & df_SOC_2024['Treatment.1'].astype(str).str.contains('Clay', case=False, na=False)]
df_SOC_2024 = df_SOC_2024[['Treatment.1', 'SOC (g/kg)']].copy()
df_SOC_2024 = df_SOC_2024.rename(columns={'Treatment.1': 'Treatment'})
df_SOC_2024['Sampling date'] = 'May_2024'
df_SOC_2024 = df_SOC_2024.set_index('Treatment')
df_SOC = pd.concat([df_SOC, df_SOC_2024])

# =========================================
#Si porewater
sheet_name = "Si_mmolL"
valid_codes = set(code_map.keys())

df_raw_Si = pd.read_excel(file_path_lysimeter_Si, sheet_name=sheet_name, header=None)

days_row_Si  = df_raw_Si.iloc[0, 1:]
dates_row_Si = df_raw_Si.iloc[2, 1:]

# --- valid days within time window ---
valid_days_Si = []
for day, date_val in zip(days_row_Si, dates_row_Si):
    try:
        d_num = float(day)
    except Exception:
        continue
    date_ts = pd.to_datetime(date_val, errors="coerce")
    if pd.notna(date_ts) and (start_date <= date_ts <= end_date):
        valid_days_Si.append(d_num)

valid_days_Si = sorted(set(valid_days_Si))

# --- build data table (row 3+) with numeric day columns ---
df_Si = df_raw_Si.iloc[3:, :].copy()
df_Si.columns = ["Sample"] + [float(x) for x in days_row_Si]
df_Si = df_Si.set_index("Sample")

# keep only valid day columns
df_Si = df_Si.loc[:, valid_days_Si]

# ensure numeric values
df_Si = df_Si.apply(pd.to_numeric, errors="coerce")

# --- keep only rows that belong to clay treatments (codes in sample name) ---
rows_keep = []
for name in df_Si.index:
    parts = str(name).strip().split("-")
    if any(p in valid_codes for p in parts):
        rows_keep.append(name)
df_Si = df_Si.loc[rows_keep]

# --- split by depth using sample name (NO def) ---
mask_12_5 = []
mask_25 = []
for name in df_Si.index:
    s = str(name)
    mask_12_5.append(bool(re.search(r"-(12\.5)\b", s)))
    mask_25.append(bool(re.search(r"-(25)\b", s)))

df_Si_12_5 = df_Si.loc[mask_12_5].copy()
df_Si_25   = df_Si.loc[mask_25].copy()

processed_data["Si_12.5cm"] = df_Si_12_5
processed_data["Si_25cm"]   = df_Si_25

# =========================
#Si in soil, HNO3 & CaCl2 2022 & 2023
path_dict_si = {
    'HNO3': file_path_HNO3,
    'CaCl2': file_path_CaCl2,
}

for ext_name, path in path_dict_si.items():
    df_Si_soil = pd.read_excel(path, sheet_name='Calc_figures')
    df_Si_soil.columns = df_Si_soil.columns.astype(str).str.strip()
    df_Si_soil = df_Si_soil[df_Si_soil['Sampling date'].astype(str).str.contains('May') & df_Si_soil['Treatment'].astype(str).str.contains('Clay', case=False, na=False)]
    df_Si_soil = df_Si_soil[['Treatment', 'Sampling date', 'Si (mg/kg)']].copy()
    df_Si_soil = df_Si_soil.groupby(['Treatment', 'Sampling date'], as_index=False).mean(numeric_only=True)
    df_Si_soil = df_Si_soil.set_index('Treatment')
    real_data_dict[f'Si_{ext_name}'] = df_Si_soil

# 2024
for sheet in path_dict_si.keys():
    df_Si_2024 = pd.read_excel(file_path_2024, sheet_name=sheet)
    df_Si_2024 = df_Si_2024[df_Si_2024['Treatment'].astype(str).str.contains('0-20') & df_Si_2024['Treatment.1'].astype(str).str.contains('Clay', case=False, na=False)]
    df_Si_2024 = df_Si_2024[['Treatment.1', 'Si (mg/kg)']].copy()
    df_Si_2024 = df_Si_2024.rename(columns={'Treatment.1': 'Treatment'})
    df_Si_2024['Sampling date'] = 'May_2024'
    df_Si_2024 = df_Si_2024.set_index('Treatment')
    real_data_dict[f'Si_{sheet}'] = pd.concat([real_data_dict[f'Si_{sheet}'], df_Si_2024])

print("Data import finished.")

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
dt = 1/(24*60) # [time resolution in days, 1 hour]
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
df_hourly = df_temp_KNMI.loc['2022':'2024']

df_temp = pd.read_excel(file_temp_and_rain,sheet_name=['2022', '2023', '2024']) 
df_temp = pd.concat(df_temp.values())
df_temp.set_index('Date', inplace=True)
df_temp = df_temp[~df_temp.index.duplicated(keep='first')].sort_index()
extended_start = pd.Timestamp(start_time) - pd.Timedelta(hours=1)
extended_end = pd.Timestamp(end_time) + pd.Timedelta(hours=1)
df_temp_subset = df_temp.loc[extended_start:extended_end]

temp_av = df_hourly['T'].mean()/10 # average temp

year_max = df_hourly['T'].max()/10
year_min = df_hourly['T'].min()/10
temp_ampl_yr = (year_max - year_min)/2 # [°C] yearly amplitude

daily_amplitudes = (df_hourly['T'].resample('D').max()/10 - df_hourly['T'].resample('D').min()/10)/2
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
temp_soil = temp_soil.resample('1min').interpolate(method='linear')
temp_soil = temp_soil.values

ET0 = pyEW.ET0(latitude,altitude,temp_air,temp_soil,temp_min,temp_max,wind,albedo,Zr,coastal,t_end,dt,day1)

# ==========================================
#vegetation
# ==========================================
T_v = t_end # [d] growth time
#In  addition, two seeds were sown in the middle of the plot to ensure eight plants per plot
#first maize harvest in 2022
df_maize1 = pd.read_excel(file_path_vegetation1,sheet_name='Massbl')
df_maize1 = df_maize1[df_maize1['Treatment'].str.contains('Clay', case=False, na=False)&~df_maize1['Treatment'].str.contains('k-feldspar', case=False, na=False)]
df_maize1 = df_maize1.groupby(['Treatment', 'Plant nr'])['Dry weight biomass (g)'].sum().reset_index()
df_maize1['Avg_Plant_Weight'] = df_maize1.groupby('Treatment')['Dry weight biomass (g)'].transform('mean')
df_maize1['k_v (g/m2)'] = df_maize1['Avg_Plant_Weight']*8/A_container# [g/m2] carrying capacity (8 plants per lysimeter)
df_maize1['v_in (g/m2)'] = 0.001*df_maize1['k_v (g/m2)']  # [g/m2] seed biomass
df_maize1 = df_maize1.groupby('Treatment', as_index=False).agg({'k_v (g/m2)': 'first',
                                                                'v_in (g/m2)': 'first'})

# first rye harvest in 2023
df_rye1 = pd.read_excel(file_path_vegetation2,sheet_name='Raw_data')
df_rye1 = df_rye1[df_rye1['Treatment'].str.contains('Clay', case=False, na=False)&~df_rye1['Treatment'].str.contains('k-feldspar', case=False, na=False)]
df_rye1 = df_rye1[['Treatment', 'Dry weight biomass (g)']].copy()
df_rye1['k_v (g/m2)'] = df_rye1['Dry weight biomass (g)']*8/A_container # [g/m2] carrying capacity (8 plants per lysimeter)
df_rye1['v_in (g/m2)'] = 0.001*df_rye1['k_v (g/m2)']  # [g/m2] seed biomass
df_rye1 = df_rye1.groupby('Treatment', as_index=False).agg({'k_v (g/m2)': 'first',
                                                            'v_in (g/m2)': 'first'}).reset_index(drop=True)

# second maize harvest in 2023
df_maize2 = pd.read_excel(file_path_vegetation3,sheet_name='Mass_bl')
df_maize2 = df_maize2[df_maize2['Treatment'].str.contains('Clay', case=False, na=False)&~df_maize2['Treatment'].str.contains('k-feldspar', case=False, na=False)]
df_maize2 = df_maize2.groupby(['Treatment', 'Treatment nr'])['Dry weight biomass (g)'].sum().reset_index()
df_maize2['Avg_Plant_Weight'] = df_maize2.groupby('Treatment')['Dry weight biomass (g)'].transform('mean')
df_maize2['k_v (g/m2)'] = df_maize2['Avg_Plant_Weight']*8/A_container# [g/m2] carrying capacity (8 plants per lysimeter)
df_maize2['v_in (g/m2)'] = 0.001*df_maize2['k_v (g/m2)']  # [g/m2] seed biomass
df_maize2 = df_maize2.groupby('Treatment', as_index=False).agg({'k_v (g/m2)': 'first',
                                                                'v_in (g/m2)': 'first'})

# second rye harvest in 2024
df_rye2 = pd.read_excel(file_path_vegetation4,sheet_name='Raw_data')
df_rye2 = df_rye2[df_rye2['Treatment'].str.contains('Clay', case=False, na=False)&~df_rye2['Treatment'].str.contains('k-feldspar', case=False, na=False)]
df_rye2 = df_rye2[['Treatment', 'Dry weight winter rye (g)']].copy()
df_rye2['k_v (g/m2)'] = df_rye2['Dry weight winter rye (g)']*8/A_container # [g/m2] carrying capacity (8 plants per lysimeter)
df_rye2['v_in (g/m2)'] = 0.001*df_rye2['k_v (g/m2)']  # [g/m2] seed biomass
df_rye2 = df_rye2.groupby('Treatment', as_index=False).agg({'k_v (g/m2)': 'first',
                                                            'v_in (g/m2)': 'first'})

RAI = 0.75 #[m2/m2] root area index
root_d = 0.3*1e-3 #[m] average root diameter
t0_v = 0 # [d] starting day of growth

# ==========================================
#save parameters
# ==========================================
seasons = [
    ('2022-05-23', '2022-09-19', df_maize1),
    ('2022-10-03', '2023-05-08', df_rye1),
    ('2023-05-19', '2023-09-14', df_maize2),
    ('2023-09-29', '2024-05-06', df_rye2),
]
n_total = len(t)
result = {grp: {} for grp in Group}

for grp in Group:
    cursor = 0
    seg_order = []
    prev_kv = 0.0  

    for idx, (start_s, end_s, df_season) in enumerate(seasons):
        # gap
        if idx > 0:
            gap_days = (pd.Timestamp(start_s) - pd.Timestamp(seasons[idx-1][1])).days - 1
            if gap_days > 0:
                gap_steps = int(gap_days / dt)
                key = f'gap_{idx}'
                result[grp][key] = {
                    'v':     np.zeros(gap_steps),
                    'k_v':   prev_kv,  
                    'start': cursor,
                    'end':   cursor + gap_steps
                }
                seg_order.append(key)
                cursor += gap_steps

        # season
        n_days  = (pd.Timestamp(end_s) - pd.Timestamp(start_s)).days + 1
        n_steps = int(n_days / dt)

        if grp in df_season['Treatment'].values:
            row = df_season[df_season['Treatment']==grp].iloc[0]
            kv  = float(row['k_v (g/m2)'])
            v_full = pyEW.veg(
                v_in=float(row['v_in (g/m2)']),
                T_v=n_days,
                k_v=kv,
                t0_v=0,
                temp_soil=temp_soil[cursor:cursor+n_steps],
                dt=dt
            )
            actual = min(len(v_full), n_steps, n_total - cursor)
            prev_kv = kv  
        else:
            kv     = 0.0
            actual = min(n_steps, n_total - cursor)
            v_full = np.zeros(actual)

        key = f'season_{idx}'
        result[grp][key] = {
            'v':     v_full[:actual],
            'k_v':   kv,
            'start': cursor,
            'end':   cursor + actual
        }
        seg_order.append(key)
        cursor += actual

    result[grp]['seg_order'] = seg_order

# ==========================================
#rain [m]
# ==========================================
rain = df_temp_subset['Rain mm'].resample('1min').interpolate(method='linear')
rain = rain.loc[start_time:end_time].values / 1000 # from mm to m, per hour and change into array

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
    day_index = int(i/(24*60))  # Assuming 60 time steps per hour
    rain[i] = rain[i] + watering[day_index]/ steps_per_day # m/day, change into m/step

# ==========================================
#moisture balance (I, Q [m], L, T, E [m/d])
# ==========================================
leach_values = [0.000569918, 0.000589896, 0.000475282, 0.000531012, 0.000478436, 0.000596205]# [m/d] leaching rate calculated from lysimeter data (leachate volume/soil area/time)
leach_by_treatment = dict(zip(Group, leach_values))

for treatment in Group:
    s_in_seg = s_in  

    for key in result[treatment]['seg_order']:
        seg = result[treatment][key]
        sl  = slice(seg['start'], seg['end'])
        n_steps_seg = seg['end'] - seg['start']

        [s, s_w, s_i, I, L, T, E, Q, Irr, n] = pyEW.moisture_balance(
            rain       = rain[sl],
            Zr         = Zr,
            soil       = soil,
            ET0        = ET0[sl],
            v          = seg['v'],
            k_v        = seg['k_v'],
            keyword_wb = keyword_wb,
            s_in       = s_in_seg,  
            t_end      = n_steps_seg,
            dt         = dt)

        seg['s']   = s
        seg['s_w'] = s_w
        seg['s_i'] = s_i
        seg['I']   = I
        seg['L']   = np.full(n_steps_seg, leach_by_treatment[treatment]*dt)
        seg['T']   = T
        seg['E']   = E
        seg['Q']   = Q
        seg['Irr'] = Irr
        seg['n']   = n

        s_in_seg = float(s[-1])  

# ==========================================
# Carbon system
# ==========================================

#Initial organic carbon
ADD = 0 # added dry litter [gOC/(m2*d)]
SOC_in = rho_bulk*13.9/1000 #[gOC/m3] initial sand OC content

#SOC balance and respiration
for treatment in Group:
    if 'Biochar' in treatment:
        Biochar_Add = rho_bulk * 5.6 / 1000
    else:
        Biochar_Add = 0

    SOC_in_seg = SOC_in  
    BIO_fast_prev = 0.0
    BIO_slow_prev = 0.0

    for idx_seg, key in enumerate(result[treatment]['seg_order']):
        seg = result[treatment][key]
        Biochar_Add = Biochar_Add if idx_seg == 0 else 0

        [SOC_total, r_het, r_aut, D, BIO_fast_prev, BIO_slow_prev,SOC_native_last] = pyEW.respiration(
            ADD=ADD, SOC_in=SOC_in_seg, CO2_air_in=CO2_air_in,
            ratio_aut_het=ratio_aut_het, soil=soil,
            s=seg['s'], v=seg['v'], k_v=seg['k_v'],
            Zr=Zr, temp_soil=temp_soil[seg['start']:seg['end']],
            dt=dt, conv_mol=conv_mol,
            Biochar_Add = Biochar_Add, k_dec_factor=0.2,
            BIO_fast_in=BIO_fast_prev,   # carry forward
            BIO_slow_in=BIO_slow_prev)   # carry forward

        seg['SOC']   = SOC_total
        seg['r_het'] = r_het
        seg['r_aut'] = r_aut
        seg['D']     = D

        SOC_in_seg = SOC_native_last


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
        'end': 22*1e-5*rho_bulk*Zr*conv_mol         
    },
    'Control - Biochar - Clay': {
        'start': CEC_tot,
        'end': 24*1e-5*rho_bulk*Zr*conv_mol           
    },
    'Olivine - Clay': {
        'start': CEC_tot,
        'end': 25*1e-5*rho_bulk*Zr*conv_mol    
    },
    'Olivine - Biochar - Clay': {
        'start': CEC_tot,
        'end': 28*1e-5*rho_bulk*Zr*conv_mol    
    },
    'Wollastonite - Clay': {
        'start': CEC_tot,
        'end': 25*1e-5*rho_bulk*Zr*conv_mol    
    },
    'Wollastonite - Biochar - Clay': {
        'start': CEC_tot,
        'end':30*1e-5*rho_bulk*Zr*conv_mol    
    }
}
print('Initial soil setting finished.')


#%% Short term calibration
# Short term calibration
# diss_f0: original values calibrated against Sept 2022 BaCl2 data (Clay-Thesis.py)
# The old *0.03 / *0.003 scaling is now replaced by the power-law decay in biogeochem.py:
#   diss_f(t) = diss_f0 * ((t_elapsed + t_ref) / (2 * t_ref)) ^ (-diss_f_n)
#   where t_ref = 119 d (Sept 2022 harvest = calibration point)
diss_f_dict = {
    'Olivine - Clay':                0.2216, #1.3%
    'Olivine - Biochar - Clay':      0.1759, #0
    'Wollastonite - Clay':           0.9353,#18.5%
    'Wollastonite - Biochar - Clay': 0.9353,#18.4%
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

#%% Calibrate diss_f_n
# Calibrate diss_f_n
# Calibrate the power-law decay exponent n for Olivine and Wollastonite
# using BaCl2 soil data at May 2023 (t~350d) and May 2024 (t~714d).
# Olivine  -> Mg signal (mg/kg)
# Wollastonite -> Ca signal (mg/kg)
# real data ：BaCl2 extraction
obs_BaCl2 = {}
for trt, ion in [('Olivine - Clay', 'Mg'),
                 ('Olivine - Biochar - Clay', 'Mg'),
                 ('Wollastonite - Clay', 'Ca'),
                 ('Wollastonite - Biochar - Clay', 'Ca')]:
    col = f'{ion} (mg/kg)'
    df_b = real_data_dict['BaCl2']
    obs_BaCl2[trt] = {
        'ion':      ion,
        'May_2023': float(df_b.loc[trt][df_b.loc[trt]['Sampling date'] == 'May_2023'][col].values[0]),
        'May_2024': float(df_b.loc[trt][df_b.loc[trt]['Sampling date'] == 'May_2024'][col].values[0]),
    }
 
# Timestep indices for May 2023 (~t=350d) and May 2024 (~t=714d)
idx_may2023 = int(350 / dt)
idx_may2024 = int(714 / dt)
 
def run_single_treatment(treatment, n_val, keyword_Si_precip=1):
    """Run full 2-year model for one treatment with a given diss_f_n value."""
    if 'Olivine' in treatment:
        mineral   = ['forsterite','enstatite','fayalite','lizardite','clinochlore']
        rock_f_in = np.array([0.54, 0.02, 0.03, 0.24, 0.06])
        SSA_in    = 3.8
        M_rock_in = 26 * rho_bulk * Zr / 1000
        psd_perc  = psd_perc_in_dunite
    elif 'Wollastonite' in treatment:
        mineral   = ['wollastonite','diopside','albite','alkali_feldspar']
        rock_f_in = np.array([0.45, 0.13, 0.09, 0.15])
        SSA_in    = 0.9
        M_rock_in = 26 * rho_bulk * Zr / 1000
        psd_perc  = psd_perc_in_wollastonite
    else:
        return None
 
    CEC_full      = linear_CEC(CEC_endpoints[treatment]['start'],
                               CEC_endpoints[treatment]['end'], n_total)
    prev_conc_in  = result[treatment]['conc_in']
    prev_f_CEC_in = result[treatment]['f_CEC_in']
    prev_d = prev_psd = prev_rock_f = prev_M_rock = prev_M_min = None
    prev_Si_in    = Si_in
    prev_CaCO3_in = CaCO3_in
    prev_MgCO3_in = MgCO3_in
    prev_pH_in    = pH_in
    all_data      = None
    first_seg     = True
    t_elapsed     = 0.0
 
    for key in result[treatment]['seg_order']:
        seg      = result[treatment][key]
        gs, ge   = seg['start'], seg['end']
        seg_days = (ge - gs) * dt
        try:
            seg_data = pyEW.biogeochem_balance(
                n=seg['n'], s=seg['s'], L=seg['L'], T=seg['T'],
                I=seg['I'], v=seg['v'], k_v=seg['k_v'],
                RAI=RAI, root_d=root_d, Zr=Zr,
                r_het=seg['r_het'], r_aut=seg['r_aut'], D=seg['D'],
                temp_soil=temp_soil[gs:ge],
                pH_in=prev_pH_in, conc_in=prev_conc_in,
                f_CEC_in=prev_f_CEC_in, K_CEC=result[treatment]['K_CEC'],
                CEC_tot=CEC_full[gs:ge], Si_in=prev_Si_in,
                CaCO3_in=prev_CaCO3_in, MgCO3_in=prev_MgCO3_in,
                M_rock_in=M_rock_in if first_seg else 0, t_app=0,
                mineral=mineral, rock_f_in=rock_f_in,
                d_in=d_in, psd_perc_in=psd_perc, SSA_in=SSA_in,
                diss_f=diss_f_dict[treatment],
                diss_f_n=n_val,
                t_elapsed_start=t_elapsed,
                dt=dt, conv_Al=conv_Al, conv_mol=conv_mol,
                keyword_add=keyword_add,
                k_prec_Si=k_Si_dict[treatment] if keyword_Si_precip == 1 else 0,
                keyword_Si_precip=keyword_Si_precip,
                d_prev=prev_d, psd_prev=prev_psd,
                rock_f_prev=prev_rock_f, M_rock_prev=prev_M_rock,
                M_min_prev=prev_M_min,
            )
        except Exception:
            return None
        first_seg  = False
        t_elapsed += seg_days
        if all_data is None:
            all_data = {k: v for k, v in seg_data.items()}
        else:
            for k, v in seg_data.items():
                if k not in all_data:
                    continue
                if isinstance(all_data[k], np.ndarray) and isinstance(v, np.ndarray) and v.ndim == 1:
                    all_data[k] = np.concatenate([all_data[k], v])

        prev_conc_in  = [float(seg_data[x][-1]) for x in ['Ca','Mg','K','Na','Al_w']]
        prev_f_CEC_in = [float(seg_data[x][-1]) for x in ['f_Ca','f_Mg','f_K','f_Na','f_Al','f_H']]
        prev_Si_in    = float(seg_data['Si'][-1])
        prev_CaCO3_in = float(seg_data['CaCO3'][-1])
        prev_MgCO3_in = float(seg_data['MgCO3'][-1])
        prev_pH_in    = float(seg_data['pH'][-1])
        if M_rock_in > 0 or prev_M_rock is not None:
            prev_d      = seg_data['d'][:, -1]
            prev_psd    = seg_data['psd'][:, -1]
            prev_rock_f = seg_data['rock_f'][:, -1]
            prev_M_rock = float(seg_data['M_rock'][-1])
            prev_M_min  = seg_data['M_min'][:, -1]
        del seg_data
        gc.collect()
    return all_data
 
 
def calibrate_n(treatment):
    """Calibrate diss_f_n for one treatment using May 2023 + May 2024 BaCl2 data."""
    info   = obs_BaCl2[treatment]
    ion    = info['ion']
    MM_ion = MM_Mg if ion == 'Mg' else MM_Ca
    obs_23 = info['May_2023']
    obs_24 = info['May_2024']
    count  = [0]
 
    def objective(params):
        n_val    = params[0]
        count[0] += 1
        all_data = run_single_treatment(treatment, n_val, keyword_Si_precip=1)
        if all_data is None:
            return 1e10
        sim_23 = float(all_data[f'{ion}_tot'][idx_may2023]) * MM_ion * 1000 / M_soil
        sim_24 = float(all_data[f'{ion}_tot'][idx_may2024]) * MM_ion * 1000 / M_soil
        # normalised RMSE over two timepoints
        err = np.sqrt(((sim_23 - obs_23)**2 + (sim_24 - obs_24)**2) / 2) / (obs_23 + 1e-10)
        print(f"  [{treatment[:22]:22s}] #{count[0]:>3} n={n_val:.3f} "
              f"sim23={sim_23:.1f} obs23={obs_23:.1f} "
              f"sim24={sim_24:.1f} obs24={obs_24:.1f} err={err:.3f}")
        return err
 
    opt = differential_evolution(
        objective, bounds=[(0.0001, 3.0)],
        maxiter=15, popsize=6, tol=0.02,
        seed=42, disp=False, workers=1, polish=True
    )
    print(f"  >>> {treatment}: best n = {opt.x[0]:.4f} (err={opt.fun:.3f})\n")
    return opt.x[0]
 
 
print("\n" + "="*60)
print("Calibrating diss_f_n (power-law decay exponent)...")
print("="*60)
for trt in ['Olivine - Clay', 'Olivine - Biochar - Clay',
            'Wollastonite - Clay', 'Wollastonite - Biochar - Clay']:
    diss_f_n_dict[trt] = calibrate_n(trt)
 
print("\n=== Final diss_f_n values ===")
for k, v in diss_f_n_dict.items():
    print(f"  {k:35s}: n = {v:.4f}")

print("Calibration complete. Using pre-calibrated values.")

#%% rock application - base model (keyword_Si_precip=0) + precip model (keyword_Si_precip=1)
# rock application
# result_base   → keyword_Si_precip=0 
# result_precip → keyword_Si_precip=1 
# ============================================================
def run_all_treatments(keyword_Si_precip, diss_f_n_override=None):
    """
    diss_f_n_override: dict {treatment: n_value} to override diss_f_n_dict,
                       used during calibration of n.
    """
    data = {}
    _diss_f_n = diss_f_n_override if diss_f_n_override is not None else diss_f_n_dict
 
    for treatment in Group:
        print(f"Starting: {treatment}")
        
        if 'Olivine' in treatment:
            mineral   = ['forsterite','enstatite','fayalite','lizardite','clinochlore']
            rock_f_in = np.array([0.54, 0.02, 0.03, 0.24, 0.06])
            SSA_in    = 3.8
            M_rock_in = 26 * rho_bulk * Zr / 1000
            psd_perc  = psd_perc_in_dunite
        elif 'Wollastonite' in treatment:
            mineral   = ['wollastonite','diopside','albite','alkali_feldspar']
            rock_f_in = np.array([0.45, 0.13, 0.09, 0.15])
            SSA_in    = 0.9
            M_rock_in = 26 * rho_bulk * Zr / 1000
            psd_perc  = psd_perc_in_wollastonite
        else:
            mineral   = ['wollastonite']
            rock_f_in = np.array([0])
            SSA_in    = 0
            M_rock_in = 0
            psd_perc  = np.array([0])
 
        CEC_full = linear_CEC(CEC_endpoints[treatment]['start'],
                              CEC_endpoints[treatment]['end'], n_total)
 
        prev_conc_in  = result[treatment]['conc_in']
        prev_f_CEC_in = result[treatment]['f_CEC_in']
        prev_d        = None
        prev_psd      = None
        prev_rock_f   = None
        prev_M_rock   = None
        prev_M_min    = None
        prev_Si_in    = Si_in
        prev_CaCO3_in = CaCO3_in
        prev_MgCO3_in = MgCO3_in
        prev_pH_in    = pH_in
 
        all_data      = None
        first_seg     = True
        t_elapsed     = 0.0  # cumulative days since mineral application
 
        for key in result[treatment]['seg_order']:
            seg = result[treatment][key]
            gs  = seg['start']
            ge  = seg['end']
            seg_days = (ge - gs) * dt  # length of this segment in days
            print(f"  {key} (steps: {ge-gs}, t_elapsed: {t_elapsed:.1f}d)")
 
            try:
                seg_data = pyEW.biogeochem_balance(
                    n         = seg['n'],
                    s         = seg['s'],
                    L         = seg['L'],
                    T         = seg['T'],
                    I         = seg['I'],
                    v         = seg['v'],
                    k_v       = seg['k_v'],
                    RAI=RAI, root_d=root_d, Zr=Zr,
                    r_het     = seg['r_het'],
                    r_aut     = seg['r_aut'],
                    D         = seg['D'],
                    temp_soil = temp_soil[gs:ge],
                    pH_in     = prev_pH_in,
                    conc_in   = prev_conc_in,
                    f_CEC_in  = prev_f_CEC_in,
                    K_CEC     = result[treatment]['K_CEC'],
                    CEC_tot   = CEC_full[gs:ge],
                    Si_in     = prev_Si_in,
                    CaCO3_in  = prev_CaCO3_in,
                    MgCO3_in  = prev_MgCO3_in,
                    M_rock_in = M_rock_in if first_seg else 0,
                    t_app     = 0,
                    mineral=mineral, rock_f_in=rock_f_in,
                    d_in=d_in, psd_perc_in=psd_perc, SSA_in=SSA_in,
                    diss_f          = diss_f_dict[treatment],
                    diss_f_n        = _diss_f_n[treatment],
                    t_elapsed_start = t_elapsed,
                    dt        = dt,
                    conv_Al=conv_Al, conv_mol=conv_mol,
                    keyword_add=keyword_add,
                    k_prec_Si = k_Si_dict[treatment] if keyword_Si_precip == 1 else 0,
                    keyword_Si_precip=keyword_Si_precip,
                    d_prev      = prev_d,
                    psd_prev    = prev_psd,
                    rock_f_prev = prev_rock_f,
                    M_rock_prev = prev_M_rock,
                    M_min_prev  = prev_M_min,
                )
            except Exception as e:
                print(f"  !! FAILED at {treatment} - {key}: {e}")
                raise
            first_seg  = False
            t_elapsed += seg_days  # advance elapsed time after each segment
 
            if all_data is None:
                all_data = {k: v for k, v in seg_data.items()}
            else:
                for k, v in seg_data.items():
                    if k not in all_data:
                        continue
                    if isinstance(all_data[k], np.ndarray) and isinstance(v, np.ndarray) and v.ndim == 1:
                        all_data[k] = np.concatenate([all_data[k], v])
                
            # carry forward
            prev_conc_in = [
                float(seg_data['Ca'][-1]),
                float(seg_data['Mg'][-1]),
                float(seg_data['K'][-1]),
                float(seg_data['Na'][-1]),
                float(seg_data['Al_w'][-1]),
            ]
            prev_f_CEC_in = [
                float(seg_data['f_Ca'][-1]),
                float(seg_data['f_Mg'][-1]),
                float(seg_data['f_K'][-1]),
                float(seg_data['f_Na'][-1]),
                float(seg_data['f_Al'][-1]),
                float(seg_data['f_H'][-1]),
            ]                
            prev_Si_in    = float(seg_data['Si'][-1])
            prev_CaCO3_in = float(seg_data['CaCO3'][-1])
            prev_MgCO3_in = float(seg_data['MgCO3'][-1])
            prev_pH_in = float(seg_data['pH'][-1])
 
            if M_rock_in > 0 or prev_M_rock is not None:
                prev_d      = seg_data['d'][:, -1]
                prev_psd    = seg_data['psd'][:, -1]
                prev_rock_f = seg_data['rock_f'][:, -1]
                prev_M_rock = float(seg_data['M_rock'][-1])
                prev_M_min  = seg_data['M_min'][:, -1]
 
            del seg_data
            gc.collect()
 
        data[treatment] = all_data
        print(f'  {treatment} finished.')
 
    return data

print('Running base model (keyword_Si_precip=0)...')
biogeochem_base  = run_all_treatments(keyword_Si_precip=0)
print('Base model finished.\n')

print('Running precip model (keyword_Si_precip=1)...')
biogeochem_precip = run_all_treatments(keyword_Si_precip=1)
print('Precip model finished.\n')

data_OS  = biogeochem_precip['Olivine - Clay']
data_OBS = biogeochem_precip['Olivine - Biochar - Clay']
data_WS  = biogeochem_precip['Wollastonite - Clay']
data_WBS = biogeochem_precip['Wollastonite - Biochar - Clay']
data_CS  = biogeochem_precip['Control - Clay']
data_CBS = biogeochem_precip['Control - Biochar - Clay']

print('All treatments finished.')

#%% Change unit - base and precip
# Change unit 
result_base   = {i: {} for i in Group}
result_precip = {i: {} for i in Group}

for i in Group:
    for version, biogeochem_dict, res in [
        ('base',   biogeochem_base,   result_base),
        ('precip', biogeochem_precip, result_precip),
    ]:
        data_current = biogeochem_dict[i]
        if data_current is None:
            print(f"  WARNING: {i} data is None, skipping")
            continue

        for ion in ions:
            res[i][f'{ion}_soil'] = pyEW.mov_avg(
                data_current[f'{ion}_tot'] * MM_conv[ion] * 1000 / M_soil, win_size)
            res[i][f'{ion}_leachate'] = pyEW.mov_avg(
                data_current[f'{ion}'] * data_current['L'] * dt, win_size)
            res[i][f'{ion} porewater con'] = pyEW.mov_avg(
                data_current[f'{ion}'] / 1000, win_size)

        res[i]['DIC']     = pyEW.mov_avg(data_current['DIC'] / 1000, win_size)
        res[i]['Alk']     = pyEW.mov_avg(data_current['Alk'] / 1000, win_size)
        res[i]['pH']      = pyEW.mov_avg(data_current['pH'], win_size)
        res[i]['Alk_CEC'] = (data_current['f_Ca'] + data_current['f_Mg'] +
                             data_current['f_Na'] + data_current['f_K'])
        res[i]['Si_soil'] = pyEW.mov_avg(
            data_current['Si_tot'] * MM_Si * 1000 / M_soil, win_size)
        res[i]['f_Ca'] = data_current['f_Ca']
        res[i]['f_Mg'] = data_current['f_Mg']
        res[i]['f_Na'] = data_current['f_Na']
        res[i]['f_K']  = data_current['f_K']
        res[i]['f_Al'] = data_current['f_Al']
        res[i]['f_H']  = data_current['f_H']

print('Unit conversion finished.')

#%% SOC  figure
# SOC  figure
# ==========================================
colors = ['#8FA5C7','#4C72B0','#EDB68A','#DD8452','#F77DE8', '#F00BD9' ]      
color = {treatment: c for treatment, c in zip(Group, colors)}
 
colors_ions = {
    'Ca': '#1f77b4', 'Mg': '#ff7f0e', 'K': '#2ca02c',
    'Na': '#d62728'
}
 
obs_days = {
    'May_2022': 0,
    'May_2023': (pd.Timestamp('2023-05-08') - pd.Timestamp(start_time)).days,
    'May_2024': (pd.Timestamp('2024-05-06') - pd.Timestamp(start_time)).days,
}
 
fig, axs = plt.subplots(1, 1, figsize=(10, 6))
 
for treatment in Group:
    SOC_full = np.concatenate([
        result[treatment][key]['SOC'] 
        for key in result[treatment]['seg_order']
    ])
    result[treatment]['SOC'] = SOC_full
 
for treatment in Group:
    c = color.get(treatment, 'black')
    ls = '-' if 'Biochar' in treatment else '--'
    axs.plot(t, result[treatment]['SOC']*1000/rho_bulk, 
             color=c, linewidth=1.5, linestyle=ls, label=treatment)
    for date_label, day in obs_days.items():
        subset = df_SOC[(df_SOC.index == treatment) & 
                        (df_SOC['Sampling date'] == date_label)]
        if not subset.empty:
            axs.scatter(day, subset['SOC (g/kg)'].values[0],
                        color=c, s=80, marker='o', zorder=5,
                        edgecolors='black', linewidths=0.5)
 
axs.set_ylabel('SOC [g kg$^{-1}$]', fontsize=12)
axs.set_xlabel('Time (days)', fontsize=12)
axs.set_title('SOC-Clay', fontsize=14, fontweight='bold', pad=12)
axs.grid(True, alpha=0.3)
 
axs.legend(loc='upper left',
           bbox_to_anchor=(1.02, 1),
           fontsize=12,
           frameon=False,
           ncol=1)
 
plt.tight_layout()
plt.savefig('SOC_figure.png', dpi=150, bbox_inches='tight')
plt.show()

#%% Base cation in soil figure
# Base cation in soil figure
# ==========================================
markers_method = {'HNO3': 'o', 'CaCl2': '^', 'BaCl2': 's'}
fig2, axes = plt.subplots(1, 3, figsize=(21,7))
ax_con, ax_woll, ax_ol = axes

#Model vs Experiment - soil concentrations
for i in Group:
    if 'Control' in i:
        ax = ax_con
        group_label = 'control'           
    elif 'Wollastonite' in i:
        ax = ax_woll
        group_label = 'wollastonite'
    elif 'Olivine' in i:
        ax = ax_ol
        group_label = 'dunite'
    
    if 'Biochar' in i:
        ls = '-'  
    else:
        ls = '--'  

    for ion in ions: 
        #Model Lines
        col_name = f"{ion} (mg/kg)"
        c = colors_ions.get(ion, 'black')
        if i in result:
            y_model = result[i][f'{ion}_soil']
            ax.plot(t, y_model, color=c, linewidth=2, alpha=0.6,linestyle=ls)
    
        for method in ['HNO3', 'CaCl2', 'BaCl2']: 
            #field exp 
            if method not in real_data_dict: continue #Cacl2 did not extract Ca
            
            df = real_data_dict[method]
            marker = markers_method[method]

            if col_name not in df.columns:continue
            if 'Biochar' in i:
                fc = c          # Fill color
                ec = 'white'    # Border color
                lw = 1          # Border thickness
            else:
                fc = 'none'     
                ec = c         
                lw = 1.5       
                
            if 'Start soil sand' in df.index:
                val_start = df.loc['Start soil sand', col_name]
                ax.scatter(0, val_start, color=c, marker=marker, 
                           s=80, facecolors=fc,edgecolors=ec,linewidths=lw, zorder=5)
            df_i = df[df.index == i]
            for _, row in df_i.iterrows():
                date_label = row.get('Sampling date', None)
                if date_label in obs_days:
                    val = row[col_name]
                    if not pd.isna(val):
                        ax.scatter(obs_days[date_label], val,
                                color=c, marker=marker, s=80,
                                facecolors=fc, edgecolors=ec,
                                linewidths=lw, zorder=5)

titles = ["Control (Clay)", "Wollastonite (Clay)", "Dunite (Clay)"]
for ax, title in zip(axes, titles):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)
axes[0].set_ylabel("Concentration (mg/kg soil)", fontsize=12)
legend_lines = [mlines.Line2D([], [], color=colors_ions[ion], linewidth=2.5, label=ion) for ion in ions]
legend_markers = [mlines.Line2D([], [], color='black', marker=markers_method[m], linestyle='None', 
                                markersize=9, markeredgecolor='white', label=m) for m in ['HNO3', 'CaCl2', 'BaCl2']]
legend_biochar = [
    mlines.Line2D([], [], color='gray', linestyle='--', linewidth=2, label='No Biochar',
                  marker='o', markersize=8, markerfacecolor='none', markeredgecolor='gray'),
                  
    mlines.Line2D([], [], color='gray', linestyle='-', linewidth=2, label='With Biochar',
                  marker='o', markersize=8, markerfacecolor='gray', markeredgecolor='white')
]
axes[1].legend(handles=legend_lines + legend_markers + legend_biochar, 
               loc='upper center', 
               bbox_to_anchor=(0.5, -0.15), 
               ncol=5, 
               frameon=False, fontsize=11)
plt.tight_layout()
plt.show()

#%% another figure of cation in soil
# another figure of cation in soil(comprising cations dissolved in the soil solution and those adsorbed onto soil colloids )
x = np.arange(len(Group))
width = 0.18
method_styles = {
    'HNO3': 'solid',      
    'CaCl2': 'dashed',    
    'BaCl2': 'dotted'     
}

method_offsets = {
    'HNO3': -0.22,
    'BaCl2': 0.0,
    'CaCl2': 0.22
}

for ion in ions:
    col_name = f"{ion} (mg/kg)"
    c = colors_ions.get(ion, 'black')

    # collect data
    all_values = []
    diff_data = {}
    
    for method in ['HNO3','BaCl2','CaCl2']:
        if method not in real_data_dict:
            continue
        df = real_data_dict[method]
        
        for i, g in enumerate(Group):
            if col_name not in df.columns:
                continue
            if g not in df.index:
                continue
            if 'Start soil sand' not in df.index:
                continue

            start_rows = df[(df.index == 'Start soil sand') | 
                ((df.index == g) & (df['Sampling date'] == 'May_2022'))]
            start = start_rows[col_name].values[0] if len(start_rows) > 0 else np.nan

            end_rows = df[(df.index == g) & (df['Sampling date'] == 'May_2024')]
            end = end_rows[col_name].values[0] if len(end_rows) > 0 else np.nan

            mid_rows = df[(df.index == g) & (df['Sampling date'] == 'May_2023')]
            mid = mid_rows[col_name].values[0] if len(mid_rows) > 0 else np.nan
            
            all_values.extend([start, end])
            diff_data[(method, g, i)] = {
                'start': start,
                'end': end,
                'diff': diff
            }
    
    # add model values
    ion_key = f'{ion}_soil'
    for i, g in enumerate(Group):
        if g in result:
            if ion_key in result[g] and len(result[g][ion_key]) > 0:
                y_model = result[g][ion_key][-1]
                all_values.append(y_model)
    
    # check if break axis is needed
    if all_values:
        p75 = np.percentile(all_values, 75)
        p95 = np.percentile(all_values, 95)
        need_break = p95 > 2.5 * p75
    else:
        need_break = False
    
    if need_break:
        # create two subplots with shared x-axis
        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(12, 8), 
                                                 sharex=True,
                                                 gridspec_kw={'height_ratios': [1, 3],
                                                             'hspace': 0.05})
        
        # set y-limits for the two axes
        bottom_max = p75 * 1.35  # slightly above the 75th percentile
        max_val = max(all_values)
        upper_values = [v for v in all_values if v > bottom_max]
        if upper_values:
            top_min = min(upper_values) * 0.75  # for the top plot, start a bit below the smallest upper value
            top_max = max(upper_values) * 1.55  # for the top plot, end a bit above the largest upper value
        else:
            top_min = p95
            top_max = max_val * 1.1
        
        ax_bottom.set_ylim(0, bottom_max)
        ax_top.set_ylim(top_min, top_max)
        
        # hide the spines between the two plots
        ax_top.spines['bottom'].set_visible(False)
        ax_bottom.spines['top'].set_visible(False)
        ax_top.xaxis.tick_top()
        ax_top.tick_params(labeltop=False)
        ax_bottom.xaxis.tick_bottom()
        
        # add diagonal lines to indicate the break
        d = 0.015  # diagonal line size in axes coordinates
        kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False, linewidth=2)
        
        # diagonal lines for the top plot
        ax_top.plot((-d, +d), (-d, +d), **kwargs)
        ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        
        # diagonal lines for the bottom plot
        kwargs.update(transform=ax_bottom.transAxes)
        ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
        ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
        
        axes = [ax_bottom, ax_top]
        
    else:
        # normal single plot
        fig, ax = plt.subplots(figsize=(12, 10))
        axes = [ax]
    
    # draw bars and markers on the appropriate axes
    for ax_current in axes:
        for (method, g, i), data in diff_data.items():
            start = data['start']
            end = data['end']
            diff = data['diff']
            
            bottom_val = min(start, end)
            height = diff

            ax_current.bar(x[i] + method_offsets[method],
                    height,
                    bottom=bottom_val,
                    width=width,
                    facecolor='none',
                    edgecolor='black',
                    linewidth=2,
                    linestyle=method_styles[method],
                    zorder=2)
            
            ax_current.scatter(x[i] + method_offsets[method],
                        start,
                        color='yellow',
                        marker='o',
                        s=50,
                        edgecolors='black',
                        linewidths=1.5,
                        zorder=3)
            
            if not np.isnan(mid):           
                ax_current.scatter(x[i] + method_offsets[method], mid,
                color='lightblue', marker='s', s=50,
                edgecolors='black', linewidths=1.5, zorder=3)

            ax_current.scatter(x[i] + method_offsets[method],
                        end,
                        color='purple',
                        marker='^',
                        s=50,
                        edgecolors='black',
                        linewidths=1.5,
                        zorder=3)
        
        # model values
        for i, g in enumerate(Group):
            if g in result:
                ion_key = f'{ion}_soil'
                if ion_key in result[g] and len(result[g][ion_key]) > 0:
                    y_model = result[g][ion_key][-1]
                    ax_current.scatter(x[i], y_model, marker='D', 
                               color=c,
                               s=120,
                               edgecolors='black',
                               linewidths=1.5, 
                               zorder=5)
        
        ax_current.grid(True, linestyle='--', alpha=0.3)
        ax_current.set_xticks(x)
    
    # set common labels, title, and legend
    if need_break:
        # set y-label on the left side of the figure
        fig.text(0.04, 0.5, 'Concentration (mg/kg)', 
                va='center', rotation='vertical', 
                fontsize=12, fontweight='bold')
        ax_top.set_title(f"{ion} - Experimental ranges vs Model (final)",
                        fontsize=14, fontweight='bold', pad=15)
        ax_bottom.set_xticklabels(Group, rotation=30, ha='right', fontsize=10)
        
        # legend
        legend_methods = [
            mlines.Line2D([], [], color='black', linewidth=2, 
                         linestyle=method_styles['CaCl2'], label='CaCl₂'),
            mlines.Line2D([], [], color='black', linewidth=2, 
                         linestyle=method_styles['BaCl2'], label='BaCl₂'),
            mlines.Line2D([], [], color='black', linewidth=2, 
                         linestyle=method_styles['HNO3'], label='HNO₃')
        ]
        
        legend_markers = [
            mlines.Line2D([], [], marker='D', linestyle='None', 
                         color=c, markersize=8, 
                         markeredgecolor='black', markeredgewidth=1.5,
                         label='Model'),
            mlines.Line2D([], [], marker='o', linestyle='None', 
                         color='yellow', markersize=6, 
                         markeredgecolor='black', markeredgewidth=0.5,
                         label='Start'),
            mlines.Line2D([], [], marker='^', linestyle='None', 
                         color='purple', markersize=6,
                         markeredgecolor='black', markeredgewidth=0.5,
                         label='End')
        ]
        
        ax_bottom.legend(handles=legend_methods + legend_markers,
                  loc='upper center',
                  bbox_to_anchor=(0.5, -0.25),
                  ncol=6,
                  frameon=True,
                  fontsize=10)
    else:
        # normal single plot
        y_margin = (max(all_values) - min(all_values)) * 0.1
        ax.set_ylim(max(0, min(all_values) - y_margin), 
                   max(all_values) + y_margin)
        
        ax.set_ylabel("Concentration (mg/kg)", fontsize=12, fontweight='bold')
        ax.set_title(f"{ion} - Experimental ranges vs Model (final)",
                    fontsize=14, fontweight='bold', pad=15)
        ax.set_xticklabels(Group, rotation=30, ha='right', fontsize=10)
        
        legend_methods = [
            mlines.Line2D([], [], color='black', linewidth=2, 
                         linestyle=method_styles['CaCl2'], label='CaCl₂'),
            mlines.Line2D([], [], color='black', linewidth=2, 
                         linestyle=method_styles['BaCl2'], label='BaCl₂'),
            mlines.Line2D([], [], color='black', linewidth=2, 
                         linestyle=method_styles['HNO3'], label='HNO₃')
        ]
        
        legend_markers = [
            mlines.Line2D([], [], marker='D', linestyle='None', 
                         color=c, markersize=8, 
                         markeredgecolor='black', markeredgewidth=1.5,
                         label='Model'),
            mlines.Line2D([], [], marker='o', linestyle='None', 
                         color='yellow', markersize=6, 
                         markeredgecolor='black', markeredgewidth=0.5,
                         label='Start'),
            mlines.Line2D([], [], marker='^', linestyle='None', 
                         color='purple', markersize=6,
                         markeredgecolor='black', markeredgewidth=0.5,
                         label='End'),
            mlines.Line2D([], [], marker='s', linestyle='None',
                         color='lightblue', markersize=6,
                         markeredgecolor='black', label='Mid (May 2023)')
        ]
        
        ax.legend(handles=legend_methods + legend_markers,
                  loc='upper center',
                  bbox_to_anchor=(0.5, -0.15),
                  ncol=6,
                  frameon=True,
                  fontsize=10)
              
    plt.tight_layout()
    plt.show()

#%% Base cation in porewater
# Base cation in porewater
# ==========================================
fig3, axes = plt.subplots(1, 3, figsize=(21,7))
ax_con, ax_woll, ax_du = axes
#Model vs Experiment - leaching concentrations
for i in Group:
    if 'Control' in i:
        ax = ax_con
        group_label = 'control'           
    elif 'Wollastonite' in i:
        ax = ax_woll
        group_label = 'wollastonite'
    elif 'Olivine' in i:
        ax = ax_du
        group_label = 'dunite'
    
    if 'Biochar' in i:
        ls = '-'  
    else:
        ls = '--'  

    for ion in ions: 
        #Model Lines
        c = colors_ions.get(ion, 'black')
        if i in result:
            y_model = result[i][f'{ion} pore water con'][1:]
            ax.plot(t[1:], y_model, color=c, linewidth=2, alpha=0.6,linestyle=ls)
    
#field exp
for ion in ions:
    df_ion = processed_data[f'{ion}_12cm']
    c = colors_ions[ion]
    for sample_name in df_ion.index:
        code_parts = sample_name.strip().split('-')
        treatment = None
        for code in code_parts:
            if code in code_map:
                treatment = code_map[code]
                break

        if 'Control' in treatment:
            ax = ax_con
        elif 'Wollastonite' in treatment:
            ax = ax_woll
        elif 'Olivine' in treatment:
            ax = ax_du
            
        if 'Biochar' in treatment:
            fc = c  
            ec = 'white'
        else:
            fc = 'none'  
            ec = c
       
        for col in df_ion.columns:
            day = float(col)
            value = df_ion.loc[sample_name, col]
            ax.scatter(day, value, color=c, marker='o', 
                        s=80, facecolors=fc, edgecolors=ec, 
                        linewidths=1.5, zorder=5, alpha=0.8)


titles = ["Control (Clay)", "Wollastonite (Clay)", "Dunite (Clay)"]
for ax, title in zip(axes, titles):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)

axes[0].set_ylabel("concentration (mmol/l)", fontsize=12)

legend_lines = [mlines.Line2D([], [], color=colors_ions[ion], linewidth=2.5, label=ion) for ion in ions]
legend_biochar = [
    mlines.Line2D([], [], color='gray', linestyle='--', linewidth=2, label='without Biochar',
                  marker='o', markersize=8, markerfacecolor='none', markeredgecolor='gray'),
    mlines.Line2D([], [], color='gray', linestyle='-', linewidth=2, label='with Biochar',
                  marker='o', markersize=8, markerfacecolor='gray', markeredgecolor='white')
]

axes[1].legend(handles=legend_lines + legend_biochar, loc='upper center', 
               bbox_to_anchor=(0.5, -0.15), ncol=6, frameon=False, fontsize=11)

plt.tight_layout()
plt.show()

#%% Alk in porewater figure
# Alk in porewater figure
# ========================================
fig_alk, axes = plt.subplots(1, 3, figsize=(21, 7))
ax_con, ax_woll, ax_ol = axes
alk_color = '#2ca02c'
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
    
    if i in result: #model
        y_model = result[i]['Alk'][1:]
        ax.plot(t[1:], y_model, color=alk_color, linewidth=2, 
                alpha=0.6, linestyle=ls)

#field exp        
df_Alk = processed_data['Alk_12cm']
for sample_name in df_Alk.index:
    code_parts = sample_name.strip().split('-')
    treatment = None
    for code in code_parts:
        if code in code_map:
            treatment = code_map[code]
            break
    
    if not treatment:
        continue
    
    if 'Control' in treatment:
        ax = ax_con
    elif 'Wollastonite' in treatment:
        ax = ax_woll
    elif 'Olivine' in treatment:
        ax = ax_ol
    
    fc = alk_color if 'Biochar' in treatment else 'none'
    ec = 'white' if 'Biochar' in treatment else alk_color
    
    for col in df_Alk.columns:
        try:
            day = float(col)
            value = df_Alk.loc[sample_name, col]
            if not np.isnan(value):
                ax.scatter(day, value, color=alk_color, marker='o',
                          s=80, facecolors=fc, edgecolors=ec,
                          linewidths=1.5, zorder=5, alpha=0.8)
        except:
            continue

titles = ["Control (Clay)", "Wollastonite (Clay)", "Olivine (Clay)"]
for ax, title in zip(axes, titles):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)

axes[0].set_ylabel("Alkalinity (mmol$_c$/L)", fontsize=12)

legend_model = [
    mlines.Line2D([], [], color=alk_color, linestyle='-',
                  linewidth=2.5, label='Model (with biochar)'),
    mlines.Line2D([], [], color=alk_color, linestyle='--',
                  linewidth=2.5, label='Model (without biochar)')
]
legend_field = [
    mlines.Line2D([], [], linestyle='None', marker='o',
                  markersize=9, markerfacecolor=alk_color,
                  markeredgecolor='white', label='Field exp. (with biochar)'),
    mlines.Line2D([], [], linestyle='None', marker='o',
                  markersize=9, markerfacecolor='none',
                  markeredgecolor=alk_color, label='Field exp. (without biochar)')
]


axes[1].legend(handles=legend_model + legend_field,
               loc='upper center', bbox_to_anchor=(0.5, -0.15),
               ncol=3, frameon=False, fontsize=11)

plt.tight_layout()
plt.show()

#%% DIC figure
# DIC 
# ==========================================
fig_dic, axes = plt.subplots(1, 3, figsize=(21, 7))
ax_con, ax_woll, ax_ol = axes

dic_color = '#1f77b4'

for i in Group:
    if 'Control' in i:
        ax = ax_con
    elif 'Wollastonite' in i:
        ax = ax_woll
    elif 'Olivine' in i:
        ax = ax_ol
    else:
        continue

    if 'Biochar' in i:
        ls = '-'
    else:
        ls = '--'

    if i in result:
        y_model = result[i]['DIC']  # [mmol/L]
        ax.plot(t, y_model, color=dic_color, linewidth=2, 
                alpha=0.6, linestyle=ls)

df_DIC = processed_data['DIC_12cm']

for sample_name in df_DIC.index:
    code_parts = sample_name.strip().split('-')
    treatment = None
    for code in code_parts:
        if code in code_map:
            treatment = code_map[code]
            break
    
    if treatment is None:
        continue

    if 'Control' in treatment:
        ax = ax_con
    elif 'Wollastonite' in treatment:
        ax = ax_woll
    elif 'Olivine' in treatment:
        ax = ax_ol
    else:
        continue
    
    if 'Biochar' in treatment:
        fc = dic_color
        ec = 'white'
    else:
        fc = 'none'
        ec = dic_color
    
    for col in df_DIC.columns:
        try:
            day = float(col)
            value = df_DIC.loc[sample_name, col]
            if not np.isnan(value):
                ax.scatter(day, value, color=dic_color, marker='o',
                          s=80, facecolors=fc, edgecolors=ec,
                          linewidths=1.5, zorder=5, alpha=0.8)
        except:
            continue

titles = ["Control (Clay)", "Wollastonite (Clay)", "Olivine (Clay)"]
for ax, title in zip(axes, titles):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim([0, max(t)])

axes[0].set_ylabel("DIC (mmol/L)", fontsize=12)

legend_model = [
    mlines.Line2D([], [], color=dic_color, linestyle='-',
                  linewidth=2.5, label='Model (with biochar)'),
    mlines.Line2D([], [], color=dic_color, linestyle='--',
                  linewidth=2.5, label='Model (without biochar)')
]
legend_field = [
    mlines.Line2D([], [], linestyle='None', marker='o',
                  markersize=9, markerfacecolor=dic_color,
                  markeredgecolor='white', label='Field exp. (with biochar)'),
    mlines.Line2D([], [], linestyle='None', marker='o',
                  markersize=9, markerfacecolor='none',
                  markeredgecolor=dic_color, label='Field exp. (without biochar)')
]

axes[1].legend(handles=legend_model + legend_field,
               loc='upper center', bbox_to_anchor=(0.5, -0.15),
               ncol=3, frameon=False, fontsize=11)

plt.tight_layout()
plt.show()


# %% pH figure
# pH figure
# ==========================================
fig_pH, axes = plt.subplots(1, 3, figsize=(21, 7))
ax_con, ax_woll, ax_ol = axes

pH_color = "#f8a407"

for i in Group:
    if 'Control' in i:
        ax = ax_con
    elif 'Wollastonite' in i:
        ax = ax_woll
    elif 'Olivine' in i:
        ax = ax_ol
    else:
        continue

    if 'Biochar' in i:
        ls = '-'
    else:
        ls = '--'

    if i in result:
        y_model = result[i]['pH']  
        ax.plot(t, y_model, color=pH_color, linewidth=2, 
                alpha=0.6, linestyle=ls)

df_pH = processed_data['pH_12cm']

for sample_name in df_pH.index:
    code_parts = sample_name.strip().split('-')
    treatment = None
    for code in code_parts:
        if code in code_map:
            treatment = code_map[code]
            break
    
    if treatment is None:
        continue

    if 'Control' in treatment:
        ax = ax_con
    elif 'Wollastonite' in treatment:
        ax = ax_woll
    elif 'Olivine' in treatment:
        ax = ax_ol
    else:
        continue
    
    if 'Biochar' in treatment:
        fc = pH_color
        ec = 'white'
    else:
        fc = 'none'
        ec = pH_color
    
    for col in df_pH.columns:
        try:
            day = float(col)
            value = df_pH.loc[sample_name, col]
            if not np.isnan(value):
                ax.scatter(day, value, color=pH_color, marker='o',
                          s=80, facecolors=fc, edgecolors=ec,
                          linewidths=1.5, zorder=5, alpha=0.8)
        except:
            continue

titles = ["Control (Clay)", "Wollastonite (Clay)", "Olivine (Clay)"]
for ax, title in zip(axes, titles):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim([0, max(t)])

axes[0].set_ylabel("pH", fontsize=12)

legend_model = [
    mlines.Line2D([], [], color=pH_color, linestyle='-',
                  linewidth=2.5, label='Model (with biochar)'),
    mlines.Line2D([], [], color=pH_color, linestyle='--',
                  linewidth=2.5, label='Model (without biochar)')
]
legend_field = [
    mlines.Line2D([], [], linestyle='None', marker='o',
                  markersize=9, markerfacecolor=pH_color,
                  markeredgecolor='white', label='Field exp. (with biochar)'),
    mlines.Line2D([], [], linestyle='None', marker='o',
                  markersize=9, markerfacecolor='none',
                  markeredgecolor=pH_color, label='Field exp. (without biochar)')
]

axes[1].legend(handles=legend_model + legend_field,
               loc='upper center', bbox_to_anchor=(0.5, -0.15),
               ncol=3, frameon=False, fontsize=11)

plt.tight_layout()
plt.show()

#%% Alk in CEC 
# Alk in CEC 
# 做敏感度分析看Alk in CEC跟CEC和BS的哪个参数关系更大
fig__alk_CEC, ax = plt.subplots(figsize=(10,10))

for idx, i in enumerate(Group): 
    c = color.get(i,'black')
    result[i]['Alk_CEC_final'] = (A_container*1e-3*CEC_tot)*result[i]['Alk_CEC'][-1]
    y_model = result[i]['Alk_CEC_final']/1000  
    ax.scatter(idx, y_model, 
                color=c, 
                s=100,
                marker = 's',
                alpha=0.6,
                zorder=3,
                edgecolors='black',
                linewidths=1.5,
                label = 'Model')
    
for idx, i in enumerate(Group):
    c = color.get(i,'black')
    if i in df_extracted.index and 'alk_CEC' in df_extracted.columns:
        value = df_extracted.loc[i, 'alk_CEC']/1000
        if not pd.isna(value):
            ax.scatter(idx, value, 
                        color=c, 
                        s=100,
                        marker='^',
                        alpha=0.6,
                        zorder=3,
                        edgecolors='black',
                        linewidths=1.5,
                        label='Measured')
   
legend = [
    mlines.Line2D([], [], marker='s', linestyle='None',
                  markerfacecolor='none', markeredgecolor='black', markersize=10,
                  label='Model'),
    mlines.Line2D([], [], marker='^', linestyle='None',
                  markerfacecolor='none', markeredgecolor='black', markersize=10,
                  label='Measured')
]
ax.legend(handles=legend, loc='upper left', framealpha=0.95, fontsize=14)
ax.set_xticks(range(len(Group)))
ax.set_xticklabels(Group, rotation=30, fontsize=10)
ax.set_ylabel('Alk in CEC (mol)', fontsize=12)
ax.set_title('Alk in CEC')
plt.tight_layout()
plt.show()

#%% Fraction of CEC
# Fraction of CEC
f_Ca_list = [data_CBS['f_Ca'], data_OBS['f_Ca'], data_WBS['f_Ca']]
f_Mg_list = [data_CBS['f_Mg'], data_OBS['f_Mg'], data_WBS['f_Mg']]
f_K_list = [data_CBS['f_K'], data_OBS['f_K'], data_WBS['f_K']]
f_Na_list = [data_CBS['f_Na'], data_OBS['f_Na'], data_WBS['f_Na']]
f_Al_list = [data_CBS['f_Al'], data_OBS['f_Al'], data_WBS['f_Al']]
f_H_list = [data_CBS['f_H'], data_OBS['f_H'], data_WBS['f_H']]
titles = ['Control - Biochar - Clay', 'Olivine - Biochar - Clay', 'Wollastonite - Biochar - Clay']
figCEC = pyEW.fig_CEC(t, f_Ca_list, f_Mg_list, f_K_list, f_Na_list, f_Al_list, f_H_list, titles)

f_Ca_list = [data_CS['f_Ca'],data_OS['f_Ca'],data_WS['f_Ca']]
f_Mg_list = [data_CS['f_Mg'],data_OS['f_Mg'],data_WS['f_Mg']]
f_K_list = [data_CS['f_K'],data_OS['f_K'],data_WS['f_K']]
f_Na_list = [data_CS['f_Na'],data_OS['f_Na'],data_WS['f_Na']]
f_Al_list = [data_CS['f_Al'],data_OS['f_Al'],data_WS['f_Al']]
f_H_list = [data_CS['f_H'],data_OS['f_H'],data_WS['f_H']]
titles = ['Control - Clay', 'Olivine - Clay', 'Wollastonite - Clay']
figCEC = pyEW.fig_CEC(t, f_Ca_list, f_Mg_list, f_K_list, f_Na_list, f_Al_list, f_H_list, titles)


#%% Carbon dioxide sequestration
# Carbon dioxide sequestration  
# ==========================================
def IC_method(data_IC,data_DIC,data_CaCO3,data_MgCO3):
    delta_DIC_L = np.trapezoid(data_DIC*1000*L,t) #- np.trapezoid(data_untr['DIC']*1000*L,t) #[mol-conv/m2]
    delta_IC_s = (data_IC[-1]-data_IC[0]) #- (data_untr['IC_tot'][-1]-data_untr['IC_tot'][0])
    delta_CaCO3 = (data_CaCO3[-1]-data_CaCO3[0]) #- (data_untr['CaCO3'][-1]-data_untr['CaCO3'][0])
    delta_MgCO3 = (data_MgCO3[-1]-data_MgCO3[0]) #- (data_untr['MgCO3'][-1]-data_untr['MgCO3'][0])
    gCO2 = 44*(delta_DIC_L + delta_IC_s + delta_CaCO3 + delta_MgCO3)/(conv_mol*rho_bulk*1e-3*Zr) # [gCO2/kg_soil]
    #gCO2 = 44*(delta_DIC_L)/(conv_mol*rho_bulk*1e-3*Zr) # [gCO2/kg_soil]
    return gCO2

def Method1_EMB(cation_treatment, cation_control, cation_mineral, mineral_dose, CDR_pot):
    if cation_mineral > 0:
        Fw = (cation_treatment - cation_control) / cation_mineral
    else:
        Fw = 0

    Fw = max(0, min(Fw, 1.0))
    CDR = (Fw * mineral_dose * CDR_pot) / M_soil  # [g CO2 kg^-1 soil]
    return Fw, CDR

# Control - Clay
Mg_CS = data_CS['Mg_tot'][-1] / conv_mol + np.trapezoid(data_CS['UP_Mg'], t) / conv_mol
Ca_CS = data_CS['Ca_tot'][-1] / conv_mol + np.trapezoid(data_CS['UP_Ca'], t) / conv_mol

# Control - Biochar - Clay
Mg_CBS = data_CBS['Mg_tot'][-1] / conv_mol + np.trapezoid(data_CBS['UP_Mg'], t) / conv_mol
Ca_CBS = data_CBS['Ca_tot'][-1] / conv_mol + np.trapezoid(data_CBS['UP_Ca'], t) / conv_mol

for i in Group:
    if 'Olivine' in i:
        primary_cation = 'Mg'
        
        base_cation_conc = {
            'Ca': 3005*1e-6/40,
            'Mg': 231670*1e-6/24.3,
            'K': 520*1e-6/39.1,
            'Na': 483*1e-6/23
        }
        
        # CDR_pot for olivine/dunite
        CDR_pot = 0.87  # [g CO2 / g weathered mineral]
        data_current = data_OBS if 'Biochar' in i else data_OS
        Mg_soil_final = data_current['Mg_tot'][-1] / conv_mol
        Mg_plant_uptake = np.trapezoid(data_current['UP_Mg'], t) / conv_mol
        Mg_treatment = Mg_soil_final + Mg_plant_uptake
        if 'Biochar' in i:
            Mg_control = Mg_CBS 
        else:
            Mg_control = Mg_CS   
        M_rock_in = 29*rho_bulk*Zr/1000  # [g mineral/m2]
        Mg_mineral = M_rock_in * base_cation_conc['Mg']
        Fw, CDR_M1 = Method1_EMB(
            cation_treatment=Mg_treatment,
            cation_control=Mg_control,
            cation_mineral=Mg_mineral,
            mineral_dose=M_rock_in,
            CDR_pot=CDR_pot
        )
        Mg_leachate = np.trapezoid(result[i]['Mg_leachate']/1e3*A_container, t)
        if Mg_mineral > 0:
            leached_fraction = Mg_leachate / (Fw * Mg_mineral) if Fw > 0 else 0
            leached_fraction = min(leached_fraction, 1.0)
        else:
            leached_fraction = 0
        Leached_CDR_M1 = CDR_M1 * leached_fraction
    elif 'Wollastonite' in i:
        primary_cation = 'Ca'
        
        base_cation_conc = {
            'Ca': 171026*1e-6/40,
            'Mg': 22312*1e-6/24.3,
            'K': 18844*1e-6/39.1,
            'Na': 15431*1e-6/23
        }
        
        # CDR_pot for wollastonite
        CDR_pot = 0.76  # [g CO2 / g weathered mineral]
        data_current = data_WBS if 'Biochar' in i else data_WS
        Ca_soil_final = data_current['Ca_tot'][-1] / conv_mol
        Ca_plant_uptake = np.trapezoid(data_current['UP_Ca'], t) / conv_mol
        Ca_treatment = Ca_soil_final + Ca_plant_uptake
        if 'Biochar' in i:
            Ca_control = Ca_CBS 
        else:
            Ca_control = Ca_CS   
        M_rock_in = 29*rho_bulk*Zr/1000  # [g mineral/m2]
        Ca_mineral = M_rock_in * base_cation_conc['Ca']
        Fw, CDR_M1 = Method1_EMB(
            cation_treatment=Ca_treatment,
            cation_control=Ca_control,
            cation_mineral=Ca_mineral,
            mineral_dose=M_rock_in,
            CDR_pot=CDR_pot
        )
        Ca_leachate = np.trapezoid(result[i]['Ca_leachate']/1e3*A_container, t)
        if Ca_mineral > 0:
            leached_fraction = Ca_leachate / (Fw * Ca_mineral) if Fw > 0 else 0
            leached_fraction = min(leached_fraction, 1.0)
        else:
            leached_fraction = 0
        Leached_CDR_M1 = CDR_M1 * leached_fraction
    
    elif 'Control' in i:
        result[i].update({'CDR': 0, 'Leached_CDR': 0, 'Fw': 0})
        continue

    result[i].update({
        'CDR': CDR_M1,
        'Leached_CDR': Leached_CDR_M1,
        'Fw': Fw
    })


#Method 3
CO2_OBS_M3 = IC_method(data_OBS['IC_tot'],data_OBS['DIC'],data_OBS['CaCO3'],data_OBS['MgCO3'])
CO2_OS_M3 = IC_method(data_OS['IC_tot'],data_OS['DIC'],data_OS['CaCO3'],data_OS['MgCO3'])
CO2_WBS_M3 = IC_method(data_WBS['IC_tot'],data_WBS['DIC'],data_WBS['CaCO3'],data_WBS['MgCO3']) 
CO2_WS_M3 = IC_method(data_WS['IC_tot'],data_WS['DIC'],data_WS['CaCO3'],data_WS['MgCO3']) 
CO2_CBS_M3 = 0.0
CO2_CS_M3 = 0.0

#Observations: from te pas chapter 6
Fw_obs_dict = {
    'Olivine - Clay': 0.190,
    'Olivine - Biochar - Clay': 0.220,
    'Wollastonite - Clay': 0.917,
    'Wollastonite - Biochar - Clay': 0.859,
}

for i, Fw_obs in Fw_obs_dict.items():
    CDR_obs = Fw_obs * M_rock_in * CDR_pot / M_soil
    result[i]['CDR_obs'] = CDR_obs
    result[i]['Fw_obs'] = Fw_obs



print("Carbon dioxide sequestration finished.")

#%% Carbon dioxide sequestration figure
# Carbon dioxide sequestration figure
# ==========================================
Method1_CDR = []
Method1_Leached_CDR = []
for i in Group:
    Method1_CDR.append(result[i].get('CDR', 0))
    Method1_Leached_CDR.append(result[i].get('Leached_CDR', 0))
Method3_total = [CO2_CBS_M3, CO2_CS_M3, CO2_OBS_M3, CO2_OS_M3, CO2_WBS_M3, CO2_WS_M3]

#Method 1 figure
fig_CDR, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(Group))
width = 0.6
ax.bar(x - width/2, Method1_CDR, width/3, label='CDR (EMB)', color='#4C72B0')
ax.bar(x + width/2, Method3_total, width/3, label='CDR (IC method)', color="#FF6F59")
ax.set_ylabel('CDR (g CO$_2$ kg$^{-1}$ soil)', fontsize=12) 
ax.set_title('Comparison of CDR by different methods during first harvest', fontsize=14) 
ax.set_xticks(x) 
ax.set_xticklabels(Group, rotation=45, ha='right') 
ax.legend() 

plt.tight_layout() 
plt.show()


# %% moisture
plt.figure(figsize=(12, 4))
plt.subplot(1, 2, 1)
plt.plot(t, rain, label='Rain')
plt.plot(t, L, label='Leaching')
plt.legend()

plt.subplot(1, 2, 2)
plt.scatter(rain, L)
plt.xlabel('Rain (m/d)')
plt.ylabel('Leaching (m/d)')
plt.show()

plt.plot(t, s)
plt.ylabel('Soil moisture')
plt.xlabel('Time (days)')
# 



# %% silicate figure — endpoint comparison (two separate figures)
# ============================================================
# 阳离子时间序列图：base model vs precip model vs 观测值
# 两条模型线几乎重叠 → 说明 Si 沉淀对阳离子影响很小
# ============================================================
colors_ions = {
    'Ca': '#1f77b4', 'Mg': '#ff7f0e', 'K': '#2ca02c', 'Na': '#d62728'
}
markers_method = {'HNO3': 'o', 'CaCl2': '^', 'BaCl2': 's'}

fig2, axes = plt.subplots(1, 3, figsize=(21, 7))
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

    for ion in ions:
        c = colors_ions.get(ion, 'black')
        col_name = f"{ion} (mg/kg)"

        # Base model line（细虚线，半透明）
        if i in result_base and f'{ion}_soil' in result_base[i]:
            ax.plot(t, result_base[i][f'{ion}_soil'],
                    color=c, linewidth=1.5, alpha=0.4,
                    linestyle=ls, label='_nolegend_')

        # Precip model line（粗实线）
        if i in result_precip and f'{ion}_soil' in result_precip[i]:
            ax.plot(t, result_precip[i][f'{ion}_soil'],
                    color=c, linewidth=2.5, alpha=0.8,
                    linestyle=ls, label='_nolegend_')

        # 观测点
        for method in ['HNO3', 'CaCl2', 'BaCl2']:
            if method not in real_data_dict:
                continue
            df = real_data_dict[method]
            if col_name not in df.columns:
                continue
            marker = markers_method[method]
            fc = c if 'Biochar' in i else 'none'
            ec = 'white' if 'Biochar' in i else c
            lw = 1 if 'Biochar' in i else 1.5

            if 'Start soil sand' in df.index:
                ax.scatter(0, df.loc['Start soil sand', col_name],
                           color=c, marker=marker, s=80,
                           facecolors=fc, edgecolors=ec,
                           linewidths=lw, zorder=5)
            if i in df.index:
                ax.scatter(t[-1], df.loc[i, col_name],
                           color=c, marker=marker, s=80,
                           facecolors=fc, edgecolors=ec,
                           linewidths=lw, zorder=5)

titles = ["Control (Clay)", "Wollastonite (Clay)", "Dunite (Clay)"]
for ax, title in zip(axes, titles):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)
axes[0].set_ylabel("Concentration (mg/kg soil)", fontsize=12)

# Legend
legend_lines   = [mlines.Line2D([], [], color=colors_ions[ion], linewidth=2.5, label=ion)
                  for ion in ions]
legend_markers = [mlines.Line2D([], [], color='black', marker=markers_method[m],
                                linestyle='None', markersize=9,
                                markeredgecolor='white', label=m)
                  for m in ['HNO3', 'CaCl2', 'BaCl2']]
legend_model   = [
    mlines.Line2D([], [], color='gray', linestyle='--', linewidth=1.5,
                  alpha=0.4, label='Base model (no Si precip)'),
    mlines.Line2D([], [], color='gray', linestyle='--', linewidth=2.5,
                  alpha=0.8, label='Precip model'),
]
legend_biochar = [
    mlines.Line2D([], [], color='gray', linestyle='--', linewidth=2,
                  label='No Biochar', marker='o', markersize=8,
                  markerfacecolor='none', markeredgecolor='gray'),
    mlines.Line2D([], [], color='gray', linestyle='-', linewidth=2,
                  label='With Biochar', marker='o', markersize=8,
                  markerfacecolor='gray', markeredgecolor='white'),
]
axes[1].legend(handles=legend_lines + legend_markers + legend_model + legend_biochar,
               loc='upper center', bbox_to_anchor=(0.5, -0.15),
               ncol=5, frameon=False, fontsize=10)
plt.tight_layout()
plt.show()


#%% Si soil figure - base vs precip (core result)
# ============================================================
# ============================================================
import matplotlib.patches as mpatches

fig_si, axes_si = plt.subplots(1, 3, figsize=(21, 7))
ax_con_si, ax_woll_si, ax_ol_si = axes_si

si_color_base   = '#4878CF'   # 蓝 - base model
si_color_precip = '#D65F5F'   # 红 - precip model
si_color_obs_c  = '#6ACC65'   # 绿 - CaCl2 观测
si_color_obs_h  = '#B47CC7'   # 紫 - HNO3 观测

for i in Group:
    if 'Control' in i:
        ax = ax_con_si
    elif 'Wollastonite' in i:
        ax = ax_woll_si
    elif 'Olivine' in i:
        ax = ax_ol_si
    else:
        continue

    ls = '-' if 'Biochar' in i else '--'
    lw = 2.0 if 'Biochar' in i else 1.5

    # Base model Si
    if i in result_base and 'Si_soil' in result_base[i]:
        ax.plot(t, result_base[i]['Si_soil'],
                color=si_color_base, linewidth=lw,
                linestyle=ls, alpha=0.8, label='_nolegend_')

    # Precip model Si
    if i in result_precip and 'Si_soil' in result_precip[i]:
        ax.plot(t, result_precip[i]['Si_soil'],
                color=si_color_precip, linewidth=lw,
                linestyle=ls, alpha=0.8, label='_nolegend_')

    # CaCl2 Si
    vals_c = processed_data['Si_CaCl2'].get(i, [])
    if len(vals_c) > 0:
        fc = si_color_obs_c if 'Biochar' in i else 'none'
        ec = 'white' if 'Biochar' in i else si_color_obs_c
        ax.scatter([t[-1]] * len(vals_c), vals_c,
                   color=si_color_obs_c, marker='o', s=80,
                   facecolors=fc, edgecolors=ec,
                   linewidths=1.5, zorder=5)

    # HNO3 Si
    vals_h = processed_data['Si_HNO3'].get(i, [])
    if len(vals_h) > 0:
        fc = si_color_obs_h if 'Biochar' in i else 'none'
        ec = 'white' if 'Biochar' in i else si_color_obs_h
        ax.scatter([t[-1]] * len(vals_h), vals_h,
                   color=si_color_obs_h, marker='s', s=80,
                   facecolors=fc, edgecolors=ec,
                   linewidths=1.5, zorder=5)

    # site-specific annotation for Si reduction percentage
    if 'Si_soil' in result_base[i] and 'Si_soil' in result_precip[i]:
        v_base   = result_base[i]['Si_soil'][-1]
        v_precip = result_precip[i]['Si_soil'][-1]
        if v_base > 0:
            imp = (v_base - v_precip) / v_base * 100
            ax.text(t[-1] * 0.7, v_base * 0.92,
                    f'Si reduction:\n{imp:.0f}%',
                    fontsize=9, color='darkred', fontweight='bold',
                    ha='center')

titles_si = ["Control (Clay)", "Wollastonite (Clay)", "Dunite (Clay)"]
for ax, title in zip(axes_si, titles_si):
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.set_xlabel("Time (days)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_ylim(bottom=0)
axes_si[0].set_ylabel("Si in soil (mg/kg)", fontsize=12)

legend_si = [
    mlines.Line2D([], [], color=si_color_base,   linewidth=2, label='Base model (no Si precip)'),
    mlines.Line2D([], [], color=si_color_precip, linewidth=2, label='Precip model'),
    mlines.Line2D([], [], color=si_color_obs_c,  linewidth=0, marker='o',
                  markersize=8, label='CaCl₂ observed'),
    mlines.Line2D([], [], color=si_color_obs_h,  linewidth=0, marker='s',
                  markersize=8, label='HNO₃ observed'),
    mlines.Line2D([], [], color='gray', linestyle='--', linewidth=2,
                  label='No Biochar'),
    mlines.Line2D([], [], color='gray', linestyle='-',  linewidth=2,
                  label='With Biochar'),
]
axes_si[1].legend(handles=legend_si,
                  loc='upper center', bbox_to_anchor=(0.5, -0.15),
                  ncol=3, frameon=False, fontsize=10)
plt.tight_layout()
plt.show()

#%%
# 测试 Clay 的 s_in
s_in_clay = 0.5

seg = result['Control - Clay'][result['Control - Clay']['seg_order'][0]]

# 重新跑moisture balance
from pyEW import moisture_balance
[s_test, s_w, s_i, I, L, T, E, Q, Irr, n] = pyEW.moisture_balance(
    rain=rain[seg['start']:seg['end']],
    Zr=Zr, soil='clay', ET0=ET0[seg['start']:seg['end']],
    v=seg['v'], k_v=seg['k_v'],
    keyword_wb=keyword_wb,
    s_in=s_in_clay,
    t_end=seg['end']-seg['start'], dt=dt)

print('s前5步:', s[:5])
print('s_h for clay: 0.47')

# 测试respiration
try:
    [SOC_total, r_het, r_aut, D, BIO_fast_prev, BIO_slow_prev, SOC_native_last] = pyEW.respiration(
        ADD=0, SOC_in=rho_bulk*13.9/1000, CO2_air_in=CO2_air_in,
        ratio_aut_het=0.5, soil='clay',
        s=s_test, v=seg['v'], k_v=seg['k_v'],
        Zr=Zr, temp_soil=temp_soil[seg['start']:seg['end']],
        dt=dt, conv_mol=conv_mol,
        Biochar_Add=0, k_dec_factor=0.2,
        BIO_fast_in=0, BIO_slow_in=0)
    print('成功! SOC前5步:', SOC_total[:5])
    print('SOC后5步:', SOC_total[-5:])
except Exception as e:
    print('错误:', e)
# %%
