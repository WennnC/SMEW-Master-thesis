import pandas as pd
import glob
import os

# ==========================================
#get all file in the folder
# ==========================================
folder_path = 'C:/Users/11734/Desktop/Thesis/Data/Weather station data' 
file_pattern = os.path.join(folder_path, "*.csv") 
all_files = glob.glob(file_pattern)

print(f"find {len(all_files)} files 。")

# ==========================================
# 2. Combine all files
# ==========================================
df_list = []

for filename in all_files:
    try:
        df = pd.read_csv(filename)         
        rain_col = next((col for col in df.columns if 'Rain' in col), None)
        temp_col = next((col for col in df.columns if 'Temperature' in col), None)
        if rain_col and temp_col:
            df_subset = df[['Date', rain_col, temp_col]].copy()
            df_subset.columns = ['Date', 'Rain mm', 'Temperature °C']
            df_list.append(df_subset)
        else:
            print(f"Columns not found in {filename}. Skipping this file.")
    except Exception as e:
        print(f"read {filename} failed: {e}")

# concatenate all dataframes in the list
if df_list:
    df_merged = pd.concat(df_list, axis=0, ignore_index=True)    
    df_merged = df_merged[df_merged['Date'].astype(str).str.contains(r'\d', na=False)]
    df_merged['Date'] = pd.to_datetime(df_merged['Date'], dayfirst=True, errors='coerce')
    if df_merged['Date'].dt.tz is not None:
        df_merged['Date'] = df_merged['Date'].dt.tz_localize(None)
    df_merged['Rain mm'] = pd.to_numeric(df_merged['Rain mm'], errors='coerce')
    df_merged['Temperature °C'] = pd.to_numeric(df_merged['Temperature °C'], errors='coerce')
    df_merged = df_merged.dropna(subset=['Rain mm', 'Temperature °C'], how='all')
    df_merged = df_merged.dropna(subset=['Date'])
    df_merged = df_merged.sort_values(by='Date')
    df_merged = df_merged.drop_duplicates(subset=['Date'], keep='first')
    print("Complteted Successfully!")

    output_file = os.path.join(folder_path, "Sinderhoeve_Cleaned.xlsx")

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        unique_years = df_merged['Date'].dt.year.unique()
        
        for year in sorted(unique_years):
            data_year = df_merged[df_merged['Date'].dt.year == year]
            
            sheet_name = str(year)
            data_year.to_excel(writer, sheet_name=sheet_name)
            print(f"  already finished sheet: {sheet_name} ")
            
    print("✅")

else:
    print("❌ No valid data to process.")