#!/usr/bin/env python3
# script to convert EDIPIC-2D snapshots to tecplot file
# Dependency: construct, numpy
# Repo: https://amatbb.amat.com/projects/EM/repos/edipic-2d/browse/Tools?at=refs%2Fheads%2Fhan%2Frelase
# Author: Han Luo <han_luo@amat.com>
# Version: 2021/06/21
# Changelog:
# - 2021-05-24: rename electron number density
# - 2021-05-25: remove StrandID
# - 2021-06-21: stepping between used snapshots
#%%
from construct import Struct
from construct import Int32sl as Int
from construct import Int32sl as Chr
from construct import Float32l as Real
from construct import Float64l as Double
from construct import Int8ul as UInt
import construct
import os
import numpy as np
import logging
import argparse
from subprocess import PIPE, Popen
from distutils.spawn import find_executable
import shlex
import concurrent.futures
from collections import OrderedDict
import glob

logging.basicConfig(format='%(levelname)s %(threadName)s : %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARN) # change this line if you don't want verbose info

def run_command(command, cwd=os.getcwd()):
    c = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE, cwd=cwd)
    o, e = c.communicate()
    o = o.decode('utf8')
    e = e.decode('utf8')
    return o, e


def parse_str(byte_list, offset = 0, dtype=Int):
    """
    Parse a null terminated string from byte_list

    Argument
    ---------
    `byte_list` byte string
    `offset`    offset of parsing

    Return
    ---------
        str
        new_offset
    """
    r = ''
    new_offset = offset
    while offset < len(byte_list):
        c, offset = parse_buffer(byte_list, dtype, offset)
        r += chr(c)
        new_offset += dtype.sizeof()
        if c == 0:
            break
    return r[:-1], new_offset


def parse_buffer(byte_list, dtype, offset=0, count=1):
    """
    Parse byte_list

    `dtype` data type. It should have a method `parse`
    `count` number of this type to parse
    `offset` start of the parsing
    `byte_list` python binary array or list like object

    """
    if dtype == str:
        return parse_str(byte_list, offset)
    else:
        if type(dtype) == construct.core.Array:
            raise TypeError('dtype should be plain data type')
        if count != 1:
            dtype = construct.Array(count, dtype)
        new_offset = offset + dtype.sizeof()
        s = dtype.parse(byte_list[offset:new_offset])
        if count != 1:
            s = list(s)
        return s, new_offset


class BinaryFile:
    '''
    Context manager with getitem
    '''
    def __init__(self, *args, **kwargs):
        if len(args) == 1:
            args = (args[0], 'rb')
        self.f = open(*args, **kwargs)
        self.size = os.path.getsize(args[0]) # size of file in byte

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __getitem__(self, index):
        if type(index) == slice:
            if index.step:
                raise ValueError(index.step)
            start = index.start or 0
            self.seek(start, os.SEEK_SET)
            if index.stop is None:
                return self.read()
            else:
                return self.read(index.stop-start)
        else:
            self.seek(index, os.SEEK_SET)
            return self.read(1)

    def read(self, *args, **kwargs):
        return self.f.read(*args, **kwargs)

    def close(self):
        return self.f.close()

    def seek(self, *args, **kwargs):
        return self.f.seek(*args, **kwargs)

    def __len__(self):
        return self.size

#%%
def read_file(file):
    with BinaryFile(file) as f:
        filesize_bytes = f.size
        offset = 0
        r_N_points_x, offset = parse_buffer(f, Real, offset)
        N_points_x = int(r_N_points_x)
        N_points_y = int(filesize_bytes / 4 / (N_points_x + 1) - 1)

        # r_dummy
        x_m, offset = parse_buffer(f, Real, offset, N_points_x)
        x_m = np.array(x_m, dtype=np.float32)

        val = np.zeros((N_points_y, N_points_x))
        y_m = np.zeros(N_points_y)
        for j in range(N_points_y):
            y_m[j], offset = parse_buffer(f, Real, offset)
            f.seek(offset, os.SEEK_SET)
            val[j] = np.fromfile(f.f, dtype=np.float32, count=N_points_x)
            offset += Real.length * N_points_x
    return val, x_m, y_m


# %%
def generate_file_pattern(nion_sp):
    '''
    Generate naming pattern

    'key' is the variable name in the final tecplot file
    'value' is the corresponding file for the variable
    '''
    file_pattern = OrderedDict({
        "Potential (V)": "F_V_2D.bin",
        "Ex (V/m)": "EX_Vm_2D.bin",
        "Ey (V/m)": "EY_Vm_2D.bin",
        "Jx (A/m<sup>2</sup>)": "JXsum_Am2_2D.bin",
        "Jy (A/m<sup>2</sup>)": "JYsum_Am2_2D.bin",
        "[ e<sup>-1</sup> ] (m<sup>-3</sup>)": "Ne_m3_2D.bin",
        "Jx [ e<sup>-1</sup> ] (A/m<sup>2</sup>)": "JXe_Am2_2D.bin",
        "Jy [ e<sup>-1</sup> ] (A/m<sup>2</sup>)": "JYe_Am2_2D.bin",
        "T<sub>e</sub>x (V)": "TXe_eV_2D.bin",
        "T<sub>e</sub>y (V)": "TYe_eV_2D.bin",
        "Vx [ e<sup>-1</sup> ] (m/s)": "VXe_ms_2D.bin",
        "Vy [ e<sup>-1</sup> ] (m/s)": "VYe_ms_2D.bin",
        "Wx [ e<sup>-1</sup> ] (V)": "WXe_eV_2D.bin",
        "Wy [ e<sup>-1</sup> ] (V)": "WYe_eV_2D.bin"
    })

    for i in range(1, 1+nion_sp):
        ion_pattern = OrderedDict({
            f"[ ION_{i} ] (m<sup>-3</sup>)": "Ni_1_m3_2D.bin",
            f"Jx [ ION_{i} ] (A/m<sup>2</sup>)": f"JXi_{i}_Am2_2D.bin",
            f"Jy [ ION_{i} ] (A/m<sup>2</sup>)": f"JYi_{i}_Am2_2D.bin",
            f"Tx [ ION_{i} ] (V)": f"TXi_{i}_eV_2D.bin",
            f"Ty [ ION_{i} ] (V)": f"TYi_{i}_eV_2D.bin",
            f"Vx [ ION_{i} ] (m/s)": f"VXi_{i}_ms_2D.bin",
            f"Vy [ ION_{i} ] (m/s)": f"VYi_{i}_ms_2D.bin",
            f"Wx [ ION_{i} ] (V)": f"WXi_{i}_eV_2D.bin",
            f"Wy [ ION_{i} ] (V)": f"WYi_{i}_eV_2D.bin",
        })
        file_pattern = OrderedDict({**file_pattern, **ion_pattern})
    return file_pattern

# %%
def load_single_snapshot(ishot, file_pattern):
    '''
    Load a single snapshot
    Arguments:
        ishot: index of the snapshot
        file_pattern: "vriable name" v.s. "file for the variable"

    Return:
        mesh: mesh used in the snapshot
        this_data: a dictionary containing all the data
    '''
    this_data = {}
    for var_name, var_file in file_pattern.items():
        file = f'_{ishot:04d}_{var_file}'
        if not os.path.isfile(file):
            logger.error(f"File [{file}] not found")
        else:
            #logger.warning(f"Load [{file}] for variable [{var_name}] at snapshot [{i}]")
            this_data[var_name], x_m, y_m = read_file(file)

    mesh = [None, None]
    mesh[0], mesh[1] = np.meshgrid(x_m, y_m)
    logger.warning(f"Finish loading snapshot {ishot:04d}")
    return mesh, this_data


def load_snapshots(snapshot_index, file_pattern):
    '''
    Load multiple snapshots with multi-threading

    Arguments:
        snapshot_index: a list of snapshot id
        file_pattern: "vriable name" v.s. "file for the variable"

    Return:
        data: a list containing the data
        mesh: mesh loaded from the last snapshot
    '''
    data = [None for i in snapshot_index]
    mesh = [None, None]
    # logger.warning(f"Start loading snapshots in parallel")
    with concurrent.futures.ProcessPoolExecutor() as executer:
        jobs = {executer.submit(load_single_snapshot, ishot, file_pattern): i for i, ishot in enumerate(snapshot_index)}
        for j in concurrent.futures.as_completed(jobs):
            i = jobs[j]
            tmp, data[i] = j.result()
            if i == len(snapshot_index)-1:
                mesh = tmp
    return mesh, data


def concatenate_files(src_file, dest_file):
    '''
    Combine "src file" to "dest_file"
    '''
    import shutil
    with open(dest_file, 'wb') as fd:
        for fname in src_file:
            with open(fname, 'rb') as fs:
                shutil.copyfileobj(fs, fd)


def write_data(fid, d, fmt='{:.5e}', max_number_per_line=10, spacing=2):
    '''
    Write data for tecplot file

    Arguments:
        fid: file handler
        d:   numpy data to write
        fmt: python formating for the data
        max_number_per_line: number of numbers in each line
        spacing: space between two number
    '''
    for i, v in enumerate(d.flatten()):
        if i%max_number_per_line != 0:
            fid.write(" "*spacing)
        fid.write(fmt.format(v))
        if i%max_number_per_line == max_number_per_line-1:
            fid.write("\n")
    fid.write("\n")


def write_zone(ii, fname, var_name, mesh, data, solution_time=0.0, strand_id=1):
    if isinstance(fname, str):
        f = open(fname, 'w')
    else:
        f = fname
    # write file header
    if ii == 0:
        f.write('VARIABLES = "X (m)", "Y (m)",' + ', '.join([f'"{i}"' for i in var_name]) + '\n')

    shape = data[var_name[0]].shape
    Imax = shape[1]
    Jmax = shape[0]
    # write zone header
    f.write(f'ZONE I={Imax} J={Jmax} DATAPACKING=BLOCK SOLUTIONTIME={solution_time}')
    if ii != 0:
        f.write(f' VARSHARELIST=([1-2]=1)')
    missing_var_index = [i for i, v in enumerate(var_name) if v not in data]
    if len(missing_var_index) > 0:
        logger.warning('Variables: ' + ",".join([f'"{var_name[i]}"' for i in missing_var_index]) + f" for snapshot {i:04d}")
        temp = ','.join([f'{i:d}' for i in missing_var_index])
        f.write(f' PASSIVEVARLIST=[{temp}]\n')
    else:
        f.write('\n')
    # write zone data
    if ii == 0:
        for d in mesh:
            write_data(f, d)
    for v in var_name:
        if v in data:
            d = data[v]
            write_data(f, d)
    f.write("\n")
    if isinstance(fname, str):
        f.close()

# %%
def export_tecplot(mesh, data, snapshot_index, snapshot_time, tec_file, parallel=True):
    '''
    Export tecplot file for the data

    Argument:
        mesh: [xcoord, ycoord]
        data: a list of data for each snapshot
        snapshot_index: index of each snapshoot
        snapshot_time: time of each snapshot
        tec_file: name of the output file
        parallel: whether to export the data with multi-processing
    '''

    var_name = list(data[0].keys())
    logging.warning('Variables: ' + ','.join([f'"{i}"' for i in var_name]))

    if parallel:
        src_file = [f'{tec_file}.{i:04d}' for i in snapshot_index]
        with concurrent.futures.ProcessPoolExecutor() as executer:
            for ii, fname in enumerate(src_file):
                executer.submit(write_zone, ii, fname, var_name, mesh, data[ii], snapshot_time[ii], snapshot_index[ii])
        concatenate_files(src_file, tec_file)
        for i in src_file:
            os.remove(i)
    else:
        with open(tec_file, 'w') as f:
            for ii in range(len(snapshot_index)):
                write_zone(ii, f, var_name, mesh, data[ii], snapshot_time[ii], snapshot_index[ii])


def get_snapshot_time():
    file = '_snapmoments.dat'
    data = np.loadtxt(file, skiprows=1, usecols=range(5)).T
    return [i for i in data[1]]


# %%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='EDIPIC-2D bin to tecplot')
    parser.add_argument('-o', '--output', default='edipic_Nd.plt', help='Output file name (default: %(default)s))')
    parser.add_argument('-is', '--istart', type=int, default=1, help='First snapshot id (default: minimum)')
    parser.add_argument('-ie', '--iend', type=int, default=-1, help='Last snapshot id (default: maximum)')
    parser.add_argument('-it', '--istep', type=int, default=1, help='Step between snapshots (default: %(default)s))')
    parser.add_argument('-b', '--binary', action='store_true', help='Export binary tecplot file')
    parser.add_argument('-v', '--version', action='version', version='2021-06-21')

    ARGS = parser.parse_args()

    all_snapshot_time = get_snapshot_time()
    istart = ARGS.istart
    if istart < 1:
        logger.warning("Miniumum available snapshot is {:d}".format(1))
        istart = 1
    iend = ARGS.iend
    if len(all_snapshot_time) < iend or iend == -1:
        if iend != -1:
            logger.warning("Maximum available snapshot is {:d}".format(len(all_snapshot_time)))
        iend = len(all_snapshot_time)
    istep = ARGS.istep
    snapshot_index = [i for i in range(istart, iend+1, istep)]
    snapshot_time = [all_snapshot_time[i-1] for i in snapshot_index]

    nion_sp = len(glob.glob(f'_{istart:04d}_Ni_*_m3_2D.bin'))
    file_pattern = generate_file_pattern(nion_sp)

    output_ascii_file = ARGS.output
    output_file = output_ascii_file
    if ARGS.binary:
        output_ascii_file = output_ascii_file + '.tmp'

    mesh, data = load_snapshots(snapshot_index, file_pattern)

    export_tecplot(mesh, data, snapshot_index, snapshot_time, output_ascii_file, parallel=True)

    if ARGS.binary:
        preplot = find_executable('preplot')
        if not preplot:
            raise ValueError('Preplot is not in PATH')
        logger.warning("Converting to tecplot binary file")
        command = 'preplot "{:s}" "{:s}"'.format(output_ascii_file, output_file)
        _o, e = run_command(command)
        if len(e) > 0:
            raise ValueError('Error: {:e}'.format(e))
        os.remove(output_ascii_file)
