#!/usr/bin/env python3
# script to convert EDIPIC-2D probes to single tecplot file
# Dependency: numpy
# Repo: https://amatbb.amat.com/projects/EM/repos/edipic-2d
# Author: Han Luo <han_luo@amat.com>
# Version: 2021/06/21
# %%
from subprocess import PIPE, Popen
from distutils.spawn import find_executable
from collections import OrderedDict
import numpy as np
import concurrent.futures
import shlex
import os
import argparse
import logging

logging.basicConfig(format='%(levelname)s %(threadName)s : %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARN) # change this line if you don't want verbose info

def run_command(command, cwd=os.getcwd()):
    c = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, cwd=cwd)
    o, e = c.communicate()
    o = o.decode('utf8')
    e = e.decode('utf8')
    return o, e


def generate_probe_file_pattern(nion_sp):
    file_pattern = OrderedDict({
        "Ex (V/m)": "dim_EX_vst.dat",
        "Ey (V/m)": "dim_EY_vst.dat",
        "Potential (V)": "dim_F_vst.dat",
        "[ e ] (m<sup>-3</sup>)": "dim_Ne_vst.dat",
    })
    for i in range(1, 1+nion_sp):
        ion_pattern = OrderedDict({
            f"[ ION_{i} ] (m<sup>-3</sup>)": f"dim_Ni_{i}_vst.dat",
        })
        file_pattern = OrderedDict({**file_pattern, **ion_pattern})
    return file_pattern


def load_probe_data(file_pattern):
    '''
    Load multiple probing data files

    Arguments:
        file_pattern: {variable_name: file_name}
    Returns:
        time: time for each step of probing
        dataOutput: a list of probing data. Each element has size num_of_timesteps x num_of_variables
    '''
    varnames = list(file_pattern.keys())
    data = [None] * len(varnames)
    data[0] = np.loadtxt(file_pattern[varnames[0]])

    ntime, nprobe = data[0].shape
    nprobe -= 1
    time = data[0][:, 0]

    with concurrent.futures.ThreadPoolExecutor() as executer:
        jobs = {executer.submit(np.loadtxt, fileName): i+1 for i,
                fileName in enumerate(list(file_pattern.values())[1:])}
        for j in concurrent.futures.as_completed(jobs):
            i = jobs[j]
            data[i] = j.result()

    # for i, fileName in enumerate(list(file_pattern.values())[1:]):
    #     data[i+1] = np.loadtxt(fileName)

    # check if we have the same number of timesteps in each file
    min_num_timestep = len(time)
    has_same_timestep = True
    time_step = {i : min_num_timestep for i in list(file_pattern.values())}
    for fileName, dd in zip(file_pattern.values(), data):
        l = len(dd[:,0])
        time_step[fileName] = l
        if l != min_num_timestep:
            has_same_timestep = False
            if l < min_num_timestep:
                min_num_timestep = l
                time = dd[:, 0]
    if not has_same_timestep:
        logging.warning(f'Probing files don\'t have the same number of time step (probably because the case is running)')
        for i, j in time_step.items():
            logging.warning(f'{i:20s}: {j} time steps')

    # make each probe as a block of data
    dataOutput = []
    for i in range(nprobe):
        d = [
            dd[:min_num_timestep, i+1]
            for dd in data]
        d = np.vstack(d).T
        dataOutput.append(d)

    return time, dataOutput


def write_zone(f, time, data, probe_info, izone, var_name):
    header = f'ZONE I={len(time)}'
    if izone > 0:
        header += ' VARSHARELIST=([1]=1)'
        header = '\n' + header
    else:
        file_header = 'VARIABLES="Time (ns)",' + \
            ', '.join([f'"{i}"' for i in var_name])
        header = file_header + '\n' + header
        data = np.vstack((time, data.T)).T

    if np.isnan(probe_info['x']):
        zone_title = f"iprobe={izone+1:d}"
    else:
        zone_title = f"iprobe={izone+1:d} x={probe_info['x']:.7e} (mm) y={probe_info['y']:.7e} (mm)"

    header += ' T="' + zone_title + '"'
    np.savetxt(f, data, fmt='%.7e', header=header, delimiter=' ', comments='')


def get_probe_locs(fname='_probelocs.dat'):
    if not os.path.isfile(fname):
        logging.error(f'File {fname} is missing, probing locations are unknown')
        return []
    probe_info = np.loadtxt(fname)
    probe = [{}] * len(probe_info)
    for i in probe_info:
        probe_num = int(i[0]) - 1
        probe[probe_num] = {
            'xi': int(i[1]),
            'yi': int(i[2]),
            'x': i[3],
            'y': i[4]
        }
    return probe


# %%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='EDIPIC probing data to tecplot file')
    parser.add_argument('-o', '--output', default='edipic_probe_1d.plt',
                        help='Output file name (default: %(default)s))')
    parser.add_argument('-b', '--binary', action='store_true', help='Export binary tecplot file')
    parser.add_argument('-v', '--version', action='version', version='2021-06-21')

    ARGS = parser.parse_args()

    import glob
    nion_sp = len(glob.glob('dim_Ni_*_vst.dat'))
    file_pattern = generate_probe_file_pattern(nion_sp)
    var_names = list(file_pattern.keys())

    # load EDIPIC data
    probe_info = get_probe_locs()
    time, data = load_probe_data(file_pattern)
    if len(probe_info) == 0:
        probe_info = [
            {
            'x': np.NAN,
            'y': np.NAN
            }
        ] * len(data)
    if len(data) != len(probe_info):
        raise ValueError("Inconsistent number of probing")
    logger.info(f"Find {len(data)} probes")

    output_ascii_file = ARGS.output
    output_file = output_ascii_file
    if ARGS.binary:
        output_ascii_file = output_ascii_file + '.tmp'

    # write Tecplot data
    if os.path.isfile(output_ascii_file):
        os.remove(output_ascii_file)
    with open(output_ascii_file, 'a') as f:
        for i, (d, p) in enumerate(zip(data, probe_info)):
            write_zone(f, time, d, p, i, var_names)

    # convert to binary
    if ARGS.binary:
        preplot = find_executable('preplot')
        if not preplot:
            raise ValueError('Preplot is not in PATH')
        logger.info("Converting to tecplot binary file")
        command = 'preplot "{:s}" "{:s}"'.format(
            output_ascii_file, output_file)
        _o, e = run_command(command)
        if len(e) > 0:
            raise ValueError('Error: {:e}'.format(e))
        os.remove(output_ascii_file)
