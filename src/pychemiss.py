#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PyChEmiss: A python emissions preprocessor for WRF-Chem regional modelling 

A preprocessor to create wrfchemi file from local emission information.
It is based on his older brother AAS4WRF.ncl 
(https://github.com/alvv1986/AAS4WRF/blob/master/AGU2017_Poster.pdf)


"""

import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe

def read_local_emiss(emission_file, sep, has_header, col_names):
    '''
    Read local emission file into a pandas DataFrame.

    Parameters
    ----------
    emission_file : str
        Path of local emission_file.
    sep : str
        Delimiter of emission_file columns.
    header : Bool
        Does the file have a header?.
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


def create_dataset_per_emiss(emiss_df, emi, lat1d, lon1d, date):
    '''
    Create a xarray dataset for an emitted species by reshaping emitted
    species column from emiss_df DataFrame based on lat, lon and date
    sizes.

    Parameters
    ----------
    emiss_df : pandas DataFrame
        created by read_emission_file.
    emi : string
        name of emitted species.
    lat1d : numpy ndarray
        latitutes of emission file.
    lon1d : numpy ndarray
        longitudes of emission file.
    date : pandas DateTimeIndex
        Hourly values from start_date and end_date.

    Returns
    -------
    ds : xarray dataset
        dataset of emitted species.

    '''
    lon, lat = np.meshgrid(lon1d, lat1d[::-1])
    ds = xr.Dataset({emi: (('Time', 'south_north', 'west_east'), 
                           emiss_df[emi].values.reshape(len(date), lat.shape[0], lon.shape[1]))},
                    coords={'Time': date,
                            'lat': (('south_north', 'west_east'), lat),
                            'lon': (('south_north', 'west_east'), lon)})
    return ds


    
def total_emiss_wrfchemi(wrfchemi, pol, molecular_mass, wrfinput):
    '''
    Calculates total emission of regridded local emissions in kTn.

    Parameters
    ----------
    wrfchmei : xarray Dataset
        Local emissions reggrided in wrfinput grid.
    pol : str
        Name of emitted pollutant.
    molecular_mass : float
        Emitted pollutant molecular mass.
    wrfinput : xarray Dataset
        WRF-Chem wrfinput file.

    Returns
    -------
    total_wrf_pol : float
        Total emitted pollutant in KTn.

    '''
    wrf_cell_area_km = wrfinput.DX *wrfinput.DY / 10**6
    wrf_pol = wrfchemi[pol]
    total_wrf_pol = (wrf_pol * wrf_cell_area_km).sum() * molecular_mass / 10**9
    return total_wrf_pol


def total_emiss_emiss_input(emiss_input, pol, molecular_mass, cell_area):
    '''
    Calculates total emission of emiss_input dataset.

    Parameters
    ----------
    emiss_input : xarray Dataset
        It contains local emitted species.
    pol : str
        Name of emitted species.
    molecular_mass : float
        Emitted species molecular mass.
    cell_area : float
        cell area of emission_file.

    Returns
    -------
    total_emiss : float
        Total emitted species in KTn.

    '''
    total_emiss = (emiss_input[pol] * cell_area).sum() * (molecular_mass / 10**9)
    return total_emiss


def nearest_method(wrfinput, emiss_input, date, cell_area):
    '''
    Perform nearest value regridding.

    Parameters
    ----------
    wrfinput : xarray Dataset
        WRF-Chem wrfinput.
    emiss_input : xarray Dataset
        It contains local emitted species.
    date : pandas DateTimeIndex
        Hourly values from start_date to end_date.
    cell_area : float
        cell_area of emission_file.

    Returns
    -------
    wrfchemi : xarray dataset
        Local emissions reggridef on wrfinput grid.

    '''
    xlat = wrfinput.XLAT.values[0, :, :]
    xlon = wrfinput.XLONG.values[0, :, :]
    grid_out = xr.Dataset(coords={'Time':date,
                                  'lat': (('south_north', 'west_east'), xlat),
                                  'lon': (('south_north', 'west_east'), xlon)})
    regridder = xe.Regridder(emiss_input, grid_out, 'nearest_s2d')
    wrfchemi = regridder(emiss_input)
    
    # Calculating CO, NO and NO2 total emiss to check conservation
    emiss_co = total_emiss_emiss_input(emiss_input, 
                                       'E_CO',
                                       28,
                                       cell_area)
    emiss_no = total_emiss_emiss_input(emiss_input, 
                                       'E_NO',
                                       30,
                                       cell_area)
    emiss_no2 = total_emiss_emiss_input(emiss_input, 
                                       'E_NO2',
                                       46,
                                       cell_area)
    
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
        
    return wrfchemi




def wrfchemi_to_netcdf(wrfchemi,wrfinput, date, emiss_names):
    '''
    Prepared wrfchemi dataset to be exported to netcdf.

    Parameters
    ----------
    wrfchemi : xarray dataset
        regriddes local emissions.
    wrfinput : xarray dataset
        WRF-Chem wrfinput.
    date : pandas DateTimeIndex
        Hourly values from start_date to end_date.
    emiss_names : list
        Names of emitted species.

    Returns
    -------
    wrfchemi : xarray dataset
        wrfchemi dataset with corrected dimensions and attributes.

    '''
    wrfchemi = (wrfchemi
                .assign_coords(emissions_zdim=0)
                .expand_dims('emissions_zdim')
                .transpose('Time', 'emissions_zdim', 'south_north', 'west_east')
                .rename({'lat':'XLAT',
                         'lon':'XLONG'}))
    
    # Adding Times variable
    wrfchemi['Times'] = xr.DataArray(date
                                     .strftime("%Y-%m-%d_%H:%M:%S")
                                     .values,
                                     dims=['Time'],
                                     coords={'Time': date.values})
    # Copying global attributes
    for key, value in wrfinput.attrs.items():
        wrfchemi.attrs[key] = value
    
    # Adding a TITLE
    wrfchemi['TITLE'] = "OUTPUT FROM PYCHEMISS PREPROCESSOR"
    
     
    # Adding attributes to emitted species variables
    for emi in emiss_names:
        wrfchemi[emi].attrs['FieldType'] = 104
        wrfchemi[emi].attrs['MemoryOrder'] = 'XYZ'
        wrfchemi[emi].attrs['description'] = 'EMISSIONS'
        wrfchemi[emi].attrs['units'] = 'mol km^2 hr^-1'
        wrfchemi[emi].attrs['stagger'] = ''
        wrfchemi[emi].attrs['coordinates'] = 'XLONG XLAT'
 
      
    return wrfchemi
    
def name_wrfchemi_file(wrfinput, dates):
    '''
    Name wrfchemi file based on dates on emission file
    (io_style_emissions = 1 or 2)

    Parameters
    ----------
    wrfinput : xarray dataset
        WRF-Chem wrfinput.
    dates : pandas DateTimeIndex
        Hourly values from start_date to end_date.

    Returns
    -------
    file_name : str
        wrfchemi name.

    '''
    dom = "".join(["d0", str(wrfinput.GRID_ID)])
    if len(dates) > 24:
        date_start = wrfinput.Times.values[0].decode("UTF-8")
        file_name = "_".join(["wrfchemi", dom, date_start])
    else:
        file_name_00z = "_".join(["wrfchemi", dom, "00z"])
        file_name_12z = "_".join(["wrfchemi", dom, "12z"])
        file_name = (file_name_00z, file_name_12z)
    return file_name


if __name__ == '__main__':
    import sys
    import yaml
    if len(sys.argv) < 2:
        print('usage: python {} src/pychemiss.yml'.format(sys.argv[0]))
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
    cell_area = config["Emissions"]["cell_area"]
    
    reggrid_method = config["Reggriding"]["method"]
    
    # Reading local emission file 
    emiss_df = read_local_emiss(emission_file, sep, has_header, 
                                col_names)
    # Getting emitted species names
    emiss_names = col_names[3:]
    
    # Emission file dimensions
    n_points = n_lon * n_lat
    lon1d = emiss_df["lon"][: n_lon].values
    lat1d = emiss_df["lat"][: n_points: n_lon].values[::-1]
    date = pd.date_range(start_date, end_date, freq='H')
    
    # Transforming text into a xarray dataset    
    # We save each emission species dataset into a dict
    DS = {}
    for emi in emiss_names:
        DS[emi] = create_dataset_per_emiss(emiss_df, emi, lat1d, lon1d, date)
        
    # Merging species datasets into one dateset
    emiss_input = xr.merge(list(DS.values()))
        
    # Reading wrfinput file
    wrfinput = xr.open_dataset(wrfinput_file)
    
    # Building wrfchemi data set
    if reggrid_method == 'nearest_s2d':
        wrfchemi = nearest_method(wrfinput, emiss_input, date, cell_area)
        
    # Building wrfchemi netcdf
    wrfchemi = wrfchemi_to_netcdf(wrfchemi, wrfinput, date, emiss_names)
   
    # Naming wrfchemi file
    output_name = name_wrfchemi_file(wrfinput, date)

    # Writting wrfchemi netcdf
    if len(date) == 24:
        wrfchemi00z = wrfchemi.isel(Time=slice(0, 12))
        wrfchemi12z = wrfchemi.isel(Time=slice(12, 24))
        wrfchemi00z.to_netcdf(output_name[0],
                               encoding={"Times":{
                                   "char_dim_name": "DateStrLen"
                                   }
                                   },
                               unlimited_dims={"Time":True})
        wrfchemi12z.to_netcdf(output_name[1],
                               encoding={"Times":{
                                   "char_dim_name": "DateStrLen"
                                   }
                                   },
                               unlimited_dims={"Time":True})

    else:
        wrfchemi.to_netcdf(output_name,
                           encoding={"Times":{
                               "char_dim_name": "DateStrLen"
                               }
                               },
                           unlimited_dims={"Time":True})
    
    
    
    print("Successful completion of PyChEmiss!")
