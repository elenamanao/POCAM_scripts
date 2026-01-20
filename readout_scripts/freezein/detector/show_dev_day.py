#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm
import os
from datetime import datetime, timezone
import numpy as np
import json


def dayplot(device, day, t_bins=0.4, noback=False, cutoff=100, showdev=False, devpos=53, string='36'):

    direct = '/data/ana/Calibration/POCAM/data/freezein/detector/'
    file_path = os.path.join(direct,day)
    files = [file for file in os.listdir(file_path) if file.endswith(".i3")]
    files.sort()
    #print(files)
    save_path = os.path.join('/data/user/leidensc/upgrade_installation_ops/freeze_in/',day)
    
    
    dev_dir = '/data/ana/Calibration/POCAM/data/freezein/device/'
    dev_path = os.path.join(dev_dir,day)
    time_file = 'run_start_end_emitter_start_end_'+device+'.json'
    with open(os.path.join(dev_path,time_file), 'r') as file:
        time_stamps = json.load(file)
    
    
    log_start = datetime.strptime(time_stamps['first_log'], "%Y-%m-%d %H:%M:%S,%f").replace(tzinfo=timezone.utc)
    log_end = datetime.strptime(time_stamps['end_log'], "%Y-%m-%d %H:%M:%S,%f").replace(tzinfo=timezone.utc)
    run_start = datetime.strptime(time_stamps['run_start'], "%Y-%m-%d %H:%M:%S,%f").replace(tzinfo=timezone.utc)
    run_end =  datetime.strptime(time_stamps['run_end'], "%Y-%m-%d %H:%M:%S,%f").replace(tzinfo=timezone.utc)
    
    #print(log_start)
    #print(log_end)
    
    #print(run_start)
    #print(run_end)
    
    
    with open(os.path.join(save_path,device), 'r') as file:
        data_dict_all = json.load(file)
    
    
    data_dict = data_dict_all[string]
    
    all_times = np.array([t for module_dict in data_dict.values() for t in module_dict.keys()])
    
    index = [i for i, x in enumerate(all_times) if float(x) > 0]
    ns_h_conv = 10**(9)
    
    t_min = float(int((log_start - run_start).total_seconds()   * 1e9))/ns_h_conv
    t_max = float(int((log_end - run_start).total_seconds()  * 1e9))/ns_h_conv
    
    
    #BIN_WIDTH_NS = 0.001*10**9/ns_h_conv  # ns # for looking into one side flash (master, slave, both)
    BIN_WIDTH_NS = t_bins#0.4*10**9/ns_h_conv      # bin size for looking at the entire run
    #BIN_WIDTH_NS = 0.01*10**9/ns_h_conv
    
    time_edges = np.arange(t_min, t_max + BIN_WIDTH_NS, BIN_WIDTH_NS)
    N_time_bins = len(time_edges) - 1
    #print(N_time_bins)
    
    
    space_list = list(data_dict.keys())
    N_space = len(space_list)
    space_vals = list(map(int, space_list))
    
    
    intensity = np.zeros((N_space, N_time_bins))
    
    
    for i, y in enumerate(space_list):
        for t_ns, charge in data_dict[y].items():
            j = int((float(t_ns)/ns_h_conv - t_min) // BIN_WIDTH_NS)
            if 0 <= j < N_time_bins:
                intensity[i, j] += charge

    if noback:
        intensity[intensity < cutoff] = np.nan
    
    
    dy = np.diff(np.array(space_vals)).mean()
    space_edges = np.concatenate(([space_vals[0] - dy/2], space_vals[:-1] + dy/2, [space_vals[-1] + dy/2]))
    
    
    fig, ax = plt.subplots(figsize = (12,9), ncols = 1, nrows =1)
    plt.pcolormesh(time_edges, space_edges, intensity, shading="flat", norm=LogNorm())
    if showdev:
        plt.hlines(devpos, xmin = t_min, xmax = t_max, linestyles='dashed', color='red', alpha = 0.5, linewidth = 2)#, size = 0.5)
    plt.gca().invert_yaxis()
    plt.xlabel("Time [s]", fontsize=15)
    plt.ylabel("DOM number", fontsize=15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    cbar = plt.colorbar()
    cbar.set_label(label="Summed charge",size=15)
    ticklabs = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(ticklabs, fontsize=14)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    #plt.xlim(104, 110)
    #plt.ylim(60,11)
    #plt.xlim(100, 110)
    #plt.xlim(95, 110)
    
    #plt.savefig('./freezein_IC86_plots/0106_pocam_str88_dev78_to_str36_no_backround.png', dpi=500, bbox_inches='tight')
    plt.show()
    
    return fig,ax