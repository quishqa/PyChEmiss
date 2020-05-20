#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aas4wrf.py: Another Assimilation System for WRF-Chem python flavored.

A mass conserving preprocessor to create wrfchemi file from 
local emission information. It is based on his older brother AAS4WRF.ncl 
(https://github.com/alvv1986/AAS4WRF/blob/master/AGU2017_Poster.pdf)


"""

import os, sys
import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe
import configparser
from cdo import Cdo

# Retrieving parameters from config file
config_file = sys.argv[1]
config = configparser.ConfigParser()
config.read(config_file)


wrfinput_file = config["Input"]["wrfinput_file"]
emission_file = config["Input"]["emission_file"]

n_lon = config.getint("Emission file", "nx")
n_lat = config.getint("Emission file", "ny")
start_date = config["Time Control"]["start_date"]
end_date = config["Time Control"]["end_date"]

reggrid_method = config['Reggriding']['method']
output_name = config["Output Style"]["output_name"]

def create_dataaaray_per_emi(emiss_df, pol, lat, lon, date):
    '''
    Create a xarray dataarray for each emission species
    by reshaping the emission dataframe (emiss_df) according the lat, lon 
    and date, add emission_zdim dimension and add attribute names to lat and
    lon dimensions.

    Parameters
    ----------
    emiss_df : pandas DataFrame
        create by read_csv emission_file
    pol : string
        name of emitted pollutant.
    lat : numpy ndarray
        latitudes of emission_file.
    lon : numpy ndarray
        longitudes of emission_file.
    date : pandas DateTimeIndex
        Hourly values from start_date to end_date.

    Returns
    -------
    emi : xr.xarray
        xarray of emitted pollutant.

    '''
    
    emi = xr.DataArray(emiss_df[pol]
                       .values
                       .reshape(len(date), len(lat), len(lon)),
                       dims=['Time', 'lat', 'lon'],
                       coords={
                           'Time': date,
                           'lat': lat[::-1],
                           'lon': lon
                       })
    # Creating new dims for wrfchemi
    emi = emi.assign_coords(emissions_zdim=0)
    emi = emi.expand_dims('emissions_zdim')

    # Ordening coordinates
    emi = emi.transpose('Time', 'emissions_zdim', 'lat', 'lon')
    emi = emi.astype('float32')

    # Giving attributes to lat and lon 
    # (this is important to calculate areas using cdo gridarea)
    emi['lat'].attrs['units'] = 'degreeN'
    emi['lon'].attrs['units'] = 'degreeE'

    return emi


def cell_bound(coord):
    '''
    Calculate the cell border coordinate from emiss_file

    Parameters
    ----------
    coord : numpy ndarray
        Latitude or longitude.

    Returns
    -------
    coord_b_final : numpy ndarray
        Latitude or longitude of cell border.

    '''
    
    dx = coord[1:] - coord[0:len(coord) - 1]
    coord_b = coord[0:len(coord) - 1] + dx/2
    coord_0 = coord[0] - dx[0]/2
    coord_f = coord[-1] + dx[-1]/2  
    coord_b_final = np.insert(coord_b, 0, coord_0)
    coord_b_final = np.append(coord_b_final, coord_f)
    return coord_b_final



def conservative_method(wrf_coords, emiss_coords, emiss_input, wrfinput):
    '''
    Perform the conservative reggridding
    adding cell bounds coordinates. Also It divides by the cell
    area to ensure conservation

    Parameters
    ----------
    wrf_coords : Dict
        It contains lat and lon wrfinput cell center.
    emiss_coords : Dict
        It contains lat and lon emiss_input cell center.
    emiss_input : xarray Dataset
        It contains local emissions pollutants xarray Dataarrays.
    wrfinput : xarray Dataset
        WRF wrfinput file.

    Returns
    -------
    wrfchemi : xarray Dataset
        Local emissions conservative reggrided in wrfinput grid.

    '''

    # Stagger grid as cell borders
    xlon_u = wrfinput.XLONG_U.values
    xlat_v = wrfinput.XLAT_V.values
    
    wrf_coords['lon_b'] = xlon_u[0, 0, :]
    wrf_coords['lat_b'] = xlat_v[0, :, 0]
    emiss_coords['lon_b'] = cell_bound(emiss_input.lon.values)
    emiss_coords['lat_b'] = cell_bound(emiss_input.lat.values)
    
    regridder = xe.Regridder(emiss_coords, wrf_coords, 
                             method = 'conservative')
    
    # We calculate the area using CDO
    emiss_input.to_netcdf("emiss_input.nc")
    cdo = Cdo()
    cdo.gridarea(input="emiss_input.nc", output="emiss_input_area.nc")
    emiss_area = xr.open_dataarray("emiss_input_area.nc")
    emiss_input_a = emiss_input/ emiss_area
    
    wrfchemi_a = regridder(emiss_input_a)
    wrf_cell_area = wrfinput.DX * wrfinput.DY
    
    wrfchemi = wrfchemi_a * wrf_cell_area     
    
    regridder.clean_weight_file()
    os.remove("emiss_input.nc")
    os.remove("emiss_input_area.nc")
    return wrfchemi



def nearest_method(wrf_coords, emiss_coords, emiss_input):
    '''
    Perform the nearest value reggriding

    Parameters
    ----------
    wrf_coords : Dict
        It contains lat and lon wrfinput cell center.
    emiss_coords : Dict
        It contains lat and lon emiss_input cell center.
    emiss_input : xarray Dataset
        It contains local emissions pollutants xarray Dataarrays.

    Returns
    -------
    wrfchemi : xarray Dataset
        Local emissions reggrided in wrfinput grid.

    '''
    
    regridder = xe.Regridder(emiss_coords, wrf_coords, method = 'nearest_s2d')
    wrfchemi = regridder(emiss_input)
    regridder.clean_weight_file()
    return wrfchemi



def write_var_attributes(wrfchemi, pol):
    '''
    Add attributes to emitted pollutant dataarrays in 
    wrfchemi dataset
    

    Parameters
    ----------
    wrfchemi : xarray Dataset
        Local emissions reggrided in wrfinput grid.
    pol : string
        Emitted pollutant.

    Returns
    -------
    None.

    '''
    ''' This function add attributes values to emitted
    species'''
    wrfchemi[pol].attrs['FieldType'] = 104
    wrfchemi[pol].attrs['MemoryOrder'] = 'XYZ'
    wrfchemi[pol].attrs['description'] = 'EMISSIONS'
    wrfchemi[pol].attrs['units'] = 'mol km^2 hr^-1'
    wrfchemi[pol].attrs['stagger'] = ''
    wrfchemi[pol].attrs['coordinates'] = 'XLONG XLAT'
    
def print_conservation(wrfchemi, emiss_input, pol):
    '''
    Print how much is conserved after regridding

    Parameters
    ----------
    wrfchemi : xarray Dataset
        Local emissions reggrided in wrfinput grid.
    emiss_input : xarray Dataset
        Local emissions from emission_file.
    pol : string
        Name of emitted pollutant.

    Returns
    -------
    None.

    '''
    diff = (wrfchemi[pol].sum() / emiss_input[pol].sum()) * 100
    print("{} is conserved: {:.2f}%".format(pol, diff.values))



# Reading local emission file CBMZ/MOSAIC mechanism
emiss_names = ['E_CO', 'E_HCHO', 'E_C2H5OH', 'E_KET', 'E_NH3', 'E_XYL', 'E_TOL',
              'E_ISO', 'E_OLI', 'E_OLT', 'E_OL2', 'E_HC8', 'E_HC5', 'E_ORA2',
              'E_ETH', 'E_ALD', 'E_CSL', 'E_SO2', 'E_HC3', 'E_NO2', 'E_NO',
              'E_CH3OH', 'E_PM25I', 'E_PM25J', 'E_SO4I', 'E_SO4J', 'E_NO3I',
              'E_NO3J', 'E_ORGI', 'E_ORGJ', 'E_ECI', 'E_ECJ', 'E_SO4C', 
              'E_NO3C', 'E_ORGC', 'E_ECC']
header_name = ["i", "lon", "lat"] + emiss_names
emiss_df = pd.read_csv(emission_file, delim_whitespace=True, names=header_name)


# Emission file dimensions
n_points = n_lon * n_lat
lon1d = emiss_df["lon"][: n_lon].values
lat1d = emiss_df["lat"][: n_points: n_lon].values[::-1]
date = pd.date_range(start_date + ' 00:00', end_date + ' 23:00', freq='H')


# Transforming text into a xarray dataset, a xarray dataset is a group
# of xarrya dataarray, a xarray dataarray is N-dimensional array with 
# labeled coordinates and dimensions  

# We save each emission species xarray in this dataset
emiss_input = xr.Dataset()

for emi in emiss_names:
    emiss_input[emi] = create_dataaaray_per_emi(emiss_df, emi, 
                                                lat1d, lon1d,
                                                date)
    
# Reading wrfinput file
wrfinput = xr.open_dataset(wrfinput_file)
xlat = wrfinput.XLAT.values
xlon = wrfinput.XLONG.values

# Creating grids used for regridding
wrf_coords = {
    'lon': xlon[0, 0, :],
    'lat': xlat[0, :, 0]
    }

emiss_coords = {
    'lon': emiss_input.lon.values,
    'lat': emiss_input.lat.values
    }


# Building wrfchemi data set
if reggrid_method == 'conservative':
    wrfchemi = conservative_method(wrf_coords, emiss_coords,
                                   emiss_input, wrfinput)
elif reggrid_method == 'nearest_s2d':
    wrfchemi = nearest_method(wrf_coords, emiss_coords, emiss_input)
    

# Building wrfchemi netcdf
for key, value in wrfinput.attrs.items():
    wrfchemi.attrs[key] = value

wrfchemi['TITLE'] = "OUTPUT FROM AAS4WRF PREPROCESSOR"

for emiss in emiss_names:
    write_var_attributes(wrfchemi, emiss)
    
wrfchemi['Times'] = xr.DataArray(
    date.strftime("%Y-%m-%d_%H:%M:%S").values,
    dims=['Time'],
    coords={'Time':date.values}
    )

wrfchemi.to_netcdf(output_name,
                   encoding={"Times":{
                       "char_dim_name": "DateStrLen"
                       }
                       },
                   unlimited_dims={"Time":True})



print_conservation(wrfchemi, emiss_input, 'E_CO')
print("Successful completion of AAS4WRF!")
