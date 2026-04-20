def CO2_pot_calc(CO2_data): #Te Pas, E. E. M. et al. Front. Clim. 5, 954064 (2023
        """
    Calculate potential CO2 sequestration based on Te Pas, E. E. M. et al. Front. Clim. 5, 954064 (2023)
    
    Parameters:
    -----------
    CO2_data : dict
        Dictionary containing:
        - rock_type: str ('olivine', 'wollastonite', 'basalt', 'anorthite', 'albite')
        - Ca_leachate, Ca_soil, Ca_initial, Ca_total: float (for Ca-based minerals)
        - Mg_leachate, Mg_soil, Mg_initial, Mg_total: float (for Mg-based minerals)
        - Na_leachate, Na_soil, Na_initial, Na_total: float (for Na-based minerals)
        - M_rock: float (mineral dose applied, in grams or kg)
        - soil_weight: float (optional, default 305.58 kg for topsoil)
    
    Returns:
    --------
    dict : {'F_weathered': float, 'CO2_P': float (g CO2/kg soil)}
    """
        rock_type = ['olivine','wollastonite','basalt', 'anorthite', 'albite']
        Balance_ions =['Ca','Mg','Na']
        if CO2_data['rock_type'] == 'wollastonite' or CO2_data['rock_type'] == 'anorthite':
            F_weathered = (CO2_data['Ca_leachate'] + CO2_data['Ca_soil']-CO2_data['Ca_initial'])/CO2_data['Ca_total']
        elif CO2_data['rock_type'] == 'olivine' or CO2_data['rock_type'] == 'basalt':
            F_weathered = (CO2_data['Mg_leachate'] + CO2_data['Mg_soil']-CO2_data['Mg_initial'])/CO2_data['Mg_total']
        elif CO2_data['rock_type'] == 'albite':
            F_weathered = (CO2_data['Na_leachate'] + CO2_data['Na_soil']-CO2_data['Na_initial'])/CO2_data['Na_total']
            if 'olivine' in CO2_data['rock_type']:
                stochiometric_factor = 1.25
            elif 'wollastonite' in CO2_data['rock_type']:
                stochiometric_factor = 0.76
            elif 'basalt' in CO2_data['rock_type']:
                stochiometric_factor = 0.38 
            elif'anorthite' in CO2_data['rock_type']:
                stochiometric_factor = 0.32
            elif 'albite' in CO2_data['rock_type']:
                stochiometric_factor = 0.17
            CO2_P = F_weathered*CO2_data['M_rock']*stochiometric_factor/305.58  # [gCO2/kg_soil] topsoil average 305.28kg
        return CO2_P

def CO2_eff_calc(CO2_data):
    CO2_E = (CO2_data['alkalinity']+CO2_data['Carbonates'])*44.0095/305.58- (CO2_data['alk_ini']+CO2_data['Carbonates_ini'])*44.0095/305.58 # [gCO2/kg_soil] topsoil average 305.28kg
    return CO2_E