import numpy as np
import pandas as pd
import time
import os
import sys
from pathlib import Path

sys.path.append(os.path.abspath('..'))
from sgp4.earth_gravity import wgs72
from sgp4 import io
from sgp4.propagation_og import sgp4 as sgp4_og
from sgp4.propagation import sgp4 as sgp4_vect

import warnings
from matplotlib import pyplot as plt
plt.style.use('seaborn-whitegrid') # Set Plot Style
plt.style.use('seaborn-poster') # Set Plot Style

SGP4_PATH = Path(__file__).parent
TLE_PATH = SGP4_PATH / "SGP4-VER.TLE"
TCPPVER_PATH = SGP4_PATH / "tcppver.out"
WHICHCONST = wgs72
COL_HEADERS = ['tsince','rx','ry','rz','vx','vy','vz']
WARNING_STR = "Warning: satnum {0} doesn't match"
SUCCESS_STR = 'success'
FAIL_STR = 'fail'

def build_sat_list(tle_path = TLE_PATH):
    with open(tle_path, 'r') as tle_file:
        tle_lines = tle_file.readlines()
    sat_list = []
    for line in tle_lines:
        if line.startswith('1'):
            line1 = line
        if line.startswith('2'):
            line2 = line[0:69]
            tstart,tend,tstep = np.array(line[69:].split(),dtype=float)
            tsince = np.arange(tstart, tend, tstep)
            # To align with tcppver.out format
            if tsince[-1] != tend:
                tsince = np.hstack((tsince,tend))
            if tsince[0] != 0.0:
                tsince = np.hstack((0,tsince))
            satrec = io.twoline2rv(line1,line2,WHICHCONST)
            sat_list.append({'satnum':satrec.satnum,
                             'satrec':satrec,
                             'tsince':tsince})
    return sat_list

def build_df_from_tcppver(tcpp_path = TCPPVER_PATH):
    tcpp_df = pd.read_csv(tcpp_path,sep='^',header=None)
    xx_mask = tcpp_df[0].str.contains('xx')
    tcpp_df['satnum'] = tcpp_df[xx_mask][0]
    tcpp_df['satnum'].ffill(inplace=True)
    tcpp_df['satnum'] = tcpp_df['satnum'].str.split(expand=True)[0]\
                                         .astype(int)
    tcpp_df = tcpp_df[~xx_mask].reset_index(drop=True)
    t_r_v = tcpp_df[0].str.split(expand=True).iloc[:,0:7].astype(float)
    t_r_v.columns = COL_HEADERS
    expected_df = pd.concat([t_r_v,tcpp_df[['satnum']]],axis=1)\
                    .reset_index(drop=True)

    return expected_df

def build_df_from_og_pysgp4(sat_list):
    og_pysgp4_df_list = []
    error_dict = {}
    for sat_dict in sat_list:
        satrec = sat_dict['satrec']
        satnum = sat_dict['satnum']
        tsince = sat_dict['tsince']
        rx = np.empty(len(tsince))
        ry = np.empty(len(tsince))
        rz = np.empty(len(tsince))
        vx = np.empty(len(tsince))
        vy = np.empty(len(tsince))
        vz = np.empty(len(tsince))
        for i,t in enumerate(tsince):
            r, v = sgp4_og(satrec,t)
            rx[i],ry[i],rz[i] = r
            vx[i],vy[i],vz[i] = v
        og_pysgp4_df = pd.DataFrame(
            data=np.array([tsince,rx,ry,rz,vx,vy,vz]).T,
            columns=COL_HEADERS
            )
        og_pysgp4_df['satnum'] = satnum
        og_pysgp4_df_list.append(og_pysgp4_df)
        error_dict[satnum] = {
            'error':satrec.error,
            'error_message':satrec.error_message,
            'method':satrec.method
            }
    og_pysgp4_df = pd.concat(og_pysgp4_df_list).reset_index(drop=True)

    return og_pysgp4_df, error_dict

def build_df_from_vect_pysgp4(sat_list):
    vect_pysgp4_df_list = []
    error_dict = {}
    for sat_dict in sat_list:
        satrec = sat_dict['satrec']
        satnum = sat_dict['satnum']
        tsince = sat_dict['tsince']
        r,v = sgp4_vect(satrec,tsince)
        nan_mask = np.isnan(r)
        if nan_mask.any():
            keep_mask = ~nan_mask[0]
            tsince = tsince[keep_mask]
            rx,ry,rz = r[:,keep_mask]
            vx,vy,vz = v[:,keep_mask]
        else:
            rx,ry,rz = r
            vx,vy,vz = v
        vect_pysgp4_df = pd.DataFrame(
            data=np.array([tsince,rx,ry,rz,vx,vy,vz]).T,
            columns=COL_HEADERS
            )
        vect_pysgp4_df['satnum'] = satnum
        vect_pysgp4_df_list.append(vect_pysgp4_df)
        error_dict[satnum] = {
            'error':satrec.error,
            'error_message':satrec.error_message,
            'method':satrec.method
            }
        
    vect_pysgp4_df = pd.concat(vect_pysgp4_df_list)\
                       .reset_index(drop=True)

    return vect_pysgp4_df, error_dict

def compare_two_sgp4_df(actual_df, actual_error_dict, expected_df):
    errors = actual_error_dict
    diff_tol = 1e-07
    for satnum,actual_group in actual_df.groupby('satnum'):
        expected_group = expected_df[expected_df['satnum'] == satnum]
        len_actual = len(actual_group)
        len_expected = len(expected_group)
        if len_actual != len_expected:
            errors[satnum]['compare'] = FAIL_STR
            print(WARNING_STR.format(satnum))
            continue
        if np.allclose(actual_group,expected_group,rtol=diff_tol):
            errors[satnum]['compare'] = SUCCESS_STR
        else:
            errors[satnum]['compare'] = FAIL_STR
            print(WARNING_STR.format(satnum))
    return errors

def compare_speed_vect_og():
    ISS_TLE = """ISS (ZARYA)
1 25544U 98067A   19340.27462131  .00000658  00000-0  19507-4 0  9998
2 25544  51.6437 228.1382 0006904   6.1312  86.6408 15.50091996201950"""
    
    tstart = 0.
    tend = 2*7*24*60.
    tstep_arr = np.arange(1.,(24.*60.)+1)[-1::-1]
    
    satname,line1,line2 = ISS_TLE.split('\n')
    vect_time = np.empty(len(tstep_arr))
    og_time = np.empty(len(tstep_arr))
    len_tsince = np.empty(len(tstep_arr))

    for i,tstep in enumerate(tstep_arr):
        tsince = np.arange(tstart,tend,tstep)
        len_tsince[i] = len(tsince)
        if tsince[-1] != tend:
            tsince = np.hstack((tsince,tend))

        start = time.time()
        satrec = io.twoline2rv(line1,line2,WHICHCONST)
        r_vect, v_vect = sgp4_vect(satrec,tsince)
        end = time.time()
        vect_time[i] = end - start

        start = time.time()
        rx_arr = np.empty(len(tsince))
        ry_arr = np.empty(len(tsince))
        rz_arr = np.empty(len(tsince))
        vx_arr = np.empty(len(tsince))
        vy_arr = np.empty(len(tsince))
        vz_arr = np.empty(len(tsince))
        for j,t in enumerate(tsince):
            r, v = sgp4_og(satrec,t)
            rx_arr[j],ry_arr[j],rz_arr[j] = r
            vx_arr[j],vy_arr[j],vz_arr[j] = v
        r = np.array([rx_arr,ry_arr,rz_arr])
        v = np.array([vx_arr,vy_arr,vz_arr])
        
        
        if not (np.allclose(r,r_vect) & np.allclose(v,v_vect)):
            raise ValueError("Values don't match")
        end = time.time()
        og_time[i] = end - start
    return (len_tsince, vect_time, og_time)