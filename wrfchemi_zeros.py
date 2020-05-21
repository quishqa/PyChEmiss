#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 19:10:18 2020

@author: quishqa
"""

import xarray as xr
import numpy as np
import pandas as pd

wrfinput = xr.open_dataset("./wrfinput_d01")
start_date = "2014-10-06"
end_date = "2014-10-12"

emiss_names = ['E_CO', 'E_HCHO', 'E_C2H5OH', 'E_KET', 'E_NH3', 'E_XYL', 
               'E_TOL', 'E_ISO', 'E_OLI', 'E_OLT', 'E_OL2', 'E_HC8', 'E_HC5', 
               'E_ORA2', 'E_ETH', 'E_ALD', 'E_CSL', 'E_SO2', 'E_HC3', 'E_NO2', 
               'E_NO', 'E_CH3OH', 'E_PM25I', 'E_PM25J', 'E_SO4I', 'E_SO4J', 
               'E_NO3I', 'E_NO3J', 'E_ORGI', 'E_ORGJ', 'E_ECI', 'E_ECJ', 
               'E_SO4C', 'E_NO3C', 'E_ORGC', 'E_ECC']

date = pd.date_range(start_date + " 00:00", end_date + " 23:00", freq="H" )


emi_zero = np.zeros((len(date), 
                     1,
                     wrfinput.south_north.shape[0], 
                     wrfinput.west_east.shape[0]))



def create_zero_emi_dataarray(emi_zero, date, wrfinput):
    emi = xr.DataArray(
        emi_zero,
        dims=['Time', 'emissions_zdim', 'south_north', 'west_east'],
        coords={
            'Time': date,
            'emissions_zdim': np.array([0]),
            'south_north': wrfinput.south_north.values,
            'west_east':wrfinput.west_east.values
            }
        )
    return(emi)
    
def write_var_attributes(wrfchemi, VAR):
    ''' This function add attributes values to emitted
    species'''
    wrfchemi[VAR].attrs['FieldType'] = 104
    wrfchemi[VAR].attrs['MemoryOrder'] = 'XYZ'
    wrfchemi[VAR].attrs['description'] = 'EMISSIONS'
    wrfchemi[VAR].attrs['units'] = 'mol km^2 hr^-1'
    wrfchemi[VAR].attrs['stagger'] = ''
    wrfchemi[VAR].attrs['coordinates'] = 'XLONG XLAT'


wrfchemi = xr.Dataset()
for emiss in emiss_names:
    wrfchemi[emiss] = create_zero_emi_dataarray(emi_zero, date, wrfinput)
    write_var_attributes(wrfchemi, emiss)
    

for key, value in wrfinput.attrs.items():
    wrfchemi.attrs[key] = value


wrfchemi['TITLE'] = 'OUTPUT FROM AAS4WRF PREPROCESSOR'

wrfchemi['Times'] = xr.DataArray(
    date.strftime("%Y-%m-%d_%H:%M:%S").values,
    dims=['Time'],
    coords={'Time':date.values}
    )


wrfchemi.to_netcdf('wrfchemi_d01_zero.nc',
                   encoding={"Times":{
                       "char_dim_name": "DateStrLen"
                       }
                       },
                   unlimited_dims={"Time":True})
