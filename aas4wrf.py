#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aas4wrf.py: Another Assimilation System for WRF-Chem python flavored.

A mass conserving preprocessor to create wrfchemi file from 
local emission information. It is based on his older brother AAS4WRF.ncl 
(https://github.com/alvv1986/AAS4WRF/blob/master/AGU2017_Poster.pdf)


"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe
from cdo import Cdo

def read_local_emiss(emission_file, sep, has_header, col_names):
    '''
    Read local emission file into a pandas DataFrame

    Parameters
    ----------
    emission_file : str
        Path of local emission_file.
    sep : str
        Delimiter of emission_file columns.
    header : Bool
        Does the file have a header.
    col_names : List
        List with emission_file columns.

    Returns
    -------
    emiss_df : pandas DataFrame
        Local emissions as Pandas DataFrame.

    '''
    if has_header:
        emiss_df = pd.read_csv(emission_file, 
                               delimiter=sep)
    else:
        emiss_df = pd.read_csv(emission_file,
                               delimiter=sep,
                               names=col_names)
    return emiss_df


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

def total_emiss_emiss_input(emiss_input, pol, molecular_mass, emiss_area):
    '''
    Calculate total emission in KTn of emiss_input domain

    Parameters
    ----------
    emiss_input : xarray Dataset
        It contains local emissions pollutants xarray Dataarrays.
    pol : str
        Name of emitted pollutant.
    molecular_mass : float
        Emitted pollutant molecular mass.
    emiss_area : xarray DataArray
        Cell grid area in m^2.

    Returns
    -------
    total_pol : float
        Total emitted pollutant in KTn.

    '''
    emiss_pol = emiss_input[pol]
    total_pol = ((emiss_pol * emiss_area / 10**6).sum() * 
                 molecular_mass / 10**9)
    return total_pol
    
def total_emiss_wrfchemi(wrfchemi, pol, molecular_mass, wrfinput):
    '''
    Calculates total emission of regridded local emissions in kTon

    Parameters
    ----------
    wrfchmei : xarray Dataset
        Local emissions conservative reggrided in wrfinput grid.
    pol : str
        Name of emitted pollutant.
    molecular_mass : float
        Emitted pollutant molecular mass.
    wrfinput : xarray Dataset
        WRF wrfinput file.

    Returns
    -------
    total_wrf_pol : TYPE
        Total emitted pollutant in KTn.

    '''
    wrf_cell_area_km = wrfinput.DX *wrfinput.DY / 10**6
    wrf_pol = wrfchemi[pol]
    total_wrf_pol = (wrf_pol * wrf_cell_area_km).sum() * molecular_mass / 10**9
    return total_wrf_pol


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
    
        
    # Conservative regridding
    wrfchemi = regridder(emiss_input) 
     
    
    # Calculating CO, NO and NO2 total emiss to check conservation
    emiss_co = total_emiss_emiss_input(emiss_input, 
                                       'E_CO',
                                       28,
                                       emiss_area)
    emiss_no = total_emiss_emiss_input(emiss_input, 
                                       'E_NO',
                                       30,
                                       emiss_area)
    emiss_no2 = total_emiss_emiss_input(emiss_input, 
                                       'E_NO2',
                                       46,
                                       emiss_area)
    
    emiss_co_wrf = total_emiss_wrfchemi(wrfchemi,
                                        'E_CO',
                                        28,
                                        wrfinput)
    emiss_no_wrf = total_emiss_wrfchemi(wrfchemi,
                                        'E_NO',
                                        30,
                                        wrfinput)
    emiss_no2_wrf = total_emiss_wrfchemi(wrfchemi,
                                         'E_NO2',
                                         46,
                                         wrfinput)
    
    print("----------INPUT TOTAL EMISSIONS----------")
    print("Total CO emission = {:.2f} kTn".format(emiss_co.values))
    print("Total NO emission = {:.2f} kTn".format(emiss_no.values))
    print("Total NO2 emission = {:.2f} kTn".format(emiss_no2.values))
    
    print("----------AFTER REGRIDDING TOTAL EMISSION----------")
    print("Total CO emission after regridding = {:.2f} kTn "
          .format(emiss_co_wrf.values))
    print("Total NO emission after regridding = {:.2f} kTn "
          .format(emiss_no_wrf.values))    
    print("Total NO2 emission after regridding = {:.2f} kTn "
          .format(emiss_no2_wrf.values))

    
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
    


if __name__ == '__main__':
    import sys
    import yaml
    if len(sys.argv) < 2:
        print('usage: python {} aasf4wrf.yml'.format(sys.argv[0]))
        sys.exit()

    config_file = sys.argv[1]
        
    # Retrieving parameters from config file
    with open(config_file) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
        
    wrfinput_file = config['Input']['wrfinput_file']
    emission_file = config['Input']['emission_file']
    
    n_lon = config['Emissions']['nx']
    n_lat = config['Emissions']['ny']
    start_date = config['Emissions']['start_date']
    end_date = config['Emissions']['end_date']
    has_header = config['Emissions']['header']
    sep = config['Emissions']['sep']
    col_names = config['Emissions']['col_names']
    
    reggrid_method = config["Reggriding"]["method"]
    
    output_name = config["Output"]["output_name"]
    


    # Reading local emission file CBMZ/MOSAIC mechanism
    emiss_df = read_local_emiss(emission_file, sep, has_header, 
                                col_names)
    # Getting emitted species names
    emiss_names = col_names[3:]
    
    
    # Emission file dimensions
    n_points = n_lon * n_lat
    lon1d = emiss_df["lon"][: n_lon].values
    lat1d = emiss_df["lat"][: n_points: n_lon].values[::-1]
    date = pd.date_range(start_date, end_date, freq='H')
    
    
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
    
    
    
    print("Successful completion of AAS4WRF!")
