import argparse
import math
import numpy as np
import os
import pandas as pd
import re
import subprocess
import sys
from datetime import datetime
from functools import partial
from multiprocessing import Pool
from numpy.core.fromnumeric import mean
from scipy import stats

import epyestim
from epyestim import covid19

sys.path.append(os.path.abspath('.'))
sys.path.append(os.path.abspath('..'))

def parse_args():
    description = "Find action points for each run of a given chain"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-act",
        "--action_plan",
        type=str,
        choices=["7-day RT", "5-day RTPlus", "7-day RTPlus", "no action", "7-day hospitalized"],
        default="5-day RTPlus"
    )
    parser.add_argument(
        "-p",
        "--downsample_prob",
        type=float,
        default=0.1
    )
    parser.add_argument(
        "-se",
        "--seed",
        type=int,
        default=1234
    )
    parser.add_argument(
        "-sw",
        "--smoothing_window",
        type=int,
        default=14
    )
    parser.add_argument(
        "-rw",
        "--r_window_size",
        type=int,
        default=14
    )
    parser.add_argument(
        "-st",
        "--series_type",
        type=str,
        choices=["SS", "CLI", "hospitalized"],
        default="SS"
    )
    parser.add_argument(
        "-csv",
        "--csv",
        type=str
    )
    parser.add_argument(
        "-suf",
        "--suffix",
        type=str,
        default=""
    )
    parser.add_argument(
        "-ser",
        "--series",
        type=int
    )
    parser.add_argument(
        "-d",
        "--delay",
        type=int,
        default=3
    )

    return parser.parse_args()

def downsampAndEstim (scen_num, downsamp_seed, traj, downsample_prob, series_type,
                      action, window, smoothing_window, r_window_size, delay):
    sc_traj = traj[traj.serial == scen_num].copy()
    sc_traj['date'] = pd.to_datetime(sc_traj['date'])
    eval_dates = pd.date_range("2020-08-13", "2021-02-14", freq="1D").tolist()

    if series_type == 'SS':
        col = "new_sym"
        sc_traj['downsample'] = downsample(sc_traj[col],
                                           prob=downsample_prob,
                                           seed=downsamp_seed)
    else:
        col = "new_adm"
        sc_traj['downsample'] = sc_traj[col]

    dt_act = None
    for dt in eval_dates:
        sc_traj_dt = sc_traj[sc_traj['date'] <= dt - pd.Timedelta(delay, 'days')]
        sc_cases = sc_traj_dt.set_index('date')['downsample']
        if estimateRtAndAction(sc_cases, action=action, window=window,
                               smoothing_window = smoothing_window,
                               r_window_size = r_window_size,
                               series_type = series_type):
            dt_act = dt
            break

    return(dt_act)
    
def consecCond(cond, window):
    cond_final = cond
    cond_shift = cond
    for i in range(window-1):
        cond_shift = cond_shift.shift()
        cond_final = cond_final & cond_shift
    return cond_final

def determineActionPoint (traj, action_plan, downsample_prob, delay,
                          downsamp_seed, smoothing_window = 14,
                          r_window_size = 14, series_type = 'SS'):
    scen_nums = traj.serial.unique()
    scen_nums.sort()
    scen_nums = scen_nums
    nscen = np.max(scen_nums % 1000)

    np.random.seed(downsamp_seed)
    downsamp_seeds = np.random.randint(10e7, size=nscen)
    eval_dates = pd.date_range("2020-08-13", "2021-02-14", freq="1D").tolist()
    actions = re.findall(r"([0-9]+)-day ([a-zA-Z]+)", action_plan)[0]

    all_dt_act = []
    for sc in scen_nums:
        print(f"Evaluating scenario number {sc}")
        sc_traj = traj[traj.serial == sc].copy()
        sc_traj['date'] = pd.to_datetime(sc_traj['date'])

        if series_type == 'SS':
            col = "new_sym"
            sc_traj['downsample'] = downsample(sc_traj[col],
                                               prob=downsample_prob,
                                               seed=downsamp_seeds[sc%1000-1])
        else:
            col = "new_adm"
            sc_traj['downsample'] = sc_traj[col]

        dt_act = None
        for dt in eval_dates:
            sc_traj_dt = sc_traj[sc_traj['date'] <= dt - pd.Timedelta(delay, 'days')]
            sc_cases = sc_traj_dt.set_index('date')['downsample']
            if estimateRtAndAction(sc_cases, action=actions[1], window=int(actions[0]),
                                   smoothing_window = smoothing_window,
                                   r_window_size = r_window_size,
                                   series_type = series_type):
                dt_act = dt
                break

        all_dt_act.append(dt_act)

    df = pd.DataFrame({'scen_num': scen_nums, 'act_pt': all_dt_act})
    return df

def estimateRtAndAction(cases, action = "RT", window = 7, return_type = "bool",
                        smoothing_window = 14, r_window_size = 14,
                        series_type = "SS"):
    # $window days Rt of more than 1 or 1.05 to trigger action
    if series_type == 'SS':
        my_continuous_distro = stats.gamma(a=5.807, scale=0.948)
        sc_distro = epyestim.discrete_distrb(my_continuous_distro)
    else:
        my_continuous_distro = stats.gamma(a=3.667, scale=3.029)
        sc_distro = epyestim.discrete_distrb(my_continuous_distro)
    #sc_distro = np.array([1])
    cond = series_type != 'SS'
    rdf = covid19.r_covid(cases, delay_distribution=sc_distro,
                          smoothing_window=smoothing_window,
                          r_window_size=r_window_size,
                          n_samples = 50,
                          auto_cutoff=series_type)
    last_nday = rdf['R_mean'][-window:]
    if action == "RT":
        cond = last_nday > 1.0
    elif action == "RTPlus":
        cond = last_nday > 1.05
    if return_type == "bool":
        return cond.all()
    elif return_type == "mean":
        return mean(last_nday)
    elif return_type == "df":
        return rdf[-window:]

def downsample(cases, prob, seed, overdisperse_multiplier=1):
    if len(cases.shape) > 1:
        sizes = cases.iloc[:,0]
    else:
        sizes = cases
    np.random.seed(seed)
    downsamp = np.random.binomial(sizes, prob)
    return downsamp

def make_traj(ser, csv):
    df = pd.read_csv(csv)
    df = df[df.serial==ser]
    return df

if __name__ == '__main__':
    args = parse_args()
    traj = make_traj(args.series, args.csv)

    if args.series_type == 'hospitalized':
        df = determineHospThreshold(flow_cfg)
        outname = f"act_pt_hosp_{args.chain}.csv"
    else:
        df = determineActionPoint(traj, args.action_plan, args.downsample_prob, args.delay,
                                  args.seed, args.smoothing_window, args.r_window_size, args.series_type)
        if args.series_type == 'SS':
            outname = f"act_pt_{args.downsample_prob}_{args.suffix}_{args.smoothing_window}_{args.r_window_size}.csv"
        else:
            outname = f"act_pt_CLI_{args.suffix}_{args.smoothing_window}_{args.r_window_size}.csv"

    df['ref_date'] = '2020-09-17'
    df['ref_date'] = pd.to_datetime(df['ref_date'])
    df['act_pt'] = pd.to_datetime(df['act_pt'])
    df['time'] = (df['act_pt'] - df['ref_date']).dt.days
    df = df.drop(['ref_date'], axis='columns')

    outpath = os.path.join("intermediate", outname)
    if os.path.exists(outpath):
        df.to_csv(os.path.join("intermediate", outname), mode='a', header=False, index=False, index_label=False, sep=' ')        
    else:
        df.to_csv(os.path.join("intermediate", outname), index=False, index_label=False, sep=' ')

