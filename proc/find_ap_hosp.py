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
        "-ht",
        "--hosp_thresh",
        type=int,
        default=600,
    )
    parser.add_argument(
        "-se",
        "--seed",
        type=int,
        default=1234
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

    return parser.parse_args()

def consecCond(cond, window):
    cond_final = cond
    cond_shift = cond
    for i in range(window-1):
        cond_shift = cond_shift.shift()
        cond_final = cond_final & cond_shift
    return cond_final

def determineHospThreshold(traj, hosp_thresh):
    scen_nums = traj.serial.unique()
    scen_nums.sort()
    nscen = np.max(scen_nums)

    all_dt_act = []
    for sc in scen_nums:
        print(f"Evaluating scenario number {sc}")
        sc_traj = traj[traj.serial == sc]
        sc_traj = sc_traj[sc_traj.time > 182]

        cond = sc_traj['HOS'] >= hosp_thresh
        cond = consecCond(cond, 3)

        sc_traj1 = sc_traj[cond]
        if not len(sc_traj1) == 0:
            dt_act = sc_traj1['date'].min() + pd.Timedelta(days=1)
        else:
            dt_act = pd.to_datetime('2017-12-23')

        all_dt_act.append(dt_act)

    df = pd.DataFrame({'scen_num': scen_nums, 'act_pt': all_dt_act})
    return df

if __name__ == '__main__':
    args = parse_args()

    traj = pd.read_csv(args.csv)
    traj['ref_date'] = '2020-02-13'
    traj['ref_date'] = pd.to_datetime(traj['ref_date'])
    traj['date'] = pd.to_datetime(traj['date'])
    traj['time'] = (traj['date'] - traj['ref_date']).dt.days
    traj = traj.drop(['ref_date'], axis='columns')

    df = determineHospThreshold(traj, args.hosp_thresh)
    df['ref_date'] = '2020-09-17'
    df['ref_date'] = pd.to_datetime(df['ref_date'])
    df['act_pt'] = pd.to_datetime(df['act_pt'])
    df['time'] = (df['act_pt'] - df['ref_date']).dt.days
    df = df.drop(['ref_date'], axis='columns')

    outname = f"act_pt_hosp_{args.suffix}.csv"
    df.to_csv(os.path.join("intermediate", outname), index=False, index_label=False, sep=' ')
