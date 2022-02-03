#! /usr/bin/env python3.6

import sys
import argparse
import configparser
import subprocess as sub
from multiprocessing import Pool

# ----------- system function ----------------

def printlog(logF, lines):
    O = open(logF, 'a')
    O.write(lines+'\n')
    print(lines)
    O.close()


def multi_run_wrapper(args):
    return subproc(*args)

def subproc(CMD, logF, log_proc):
    call = 0
    if logF == "stdout":
        printlog(log_proc, f'$ {CMD}')
        call = sub.call(CMD,  shell=True,universal_newlines=True)
        pass
    else:
        printlog(log_proc, f'$ {CMD}') ## FORTEST
        call = sub.call(f'{CMD} 2> {logF}', shell=True, universal_newlines=True)
        pass
    if call == 0 :
        return 0
    else :
        printlog(log_proc, "\n*** There's a problem executing the upper command.")
        printlog(log_proc, f"*** check logfile > {logF}\n")
        return 1


def set_job_number(lib_num, given_core):
    if lib_num >= given_core:
        core =1
        job_number = given_core

    elif lib_num< given_core:
        job_number = lib_num
        core = given_core//job_number

    return job_number, core


class Parameters():
    def __init__(self, config,key1):
        for key2 in config[key1]:
            self.__dict__[key2] = config[key1][key2]

    def param_keys(self):
        return list(self.__dict__.keys())

    def param_values(self):
        return list(self.__dict__.values())


# -------   data preparence -----------------

def UCSC_download(ucsc_name, f , refD):
    global log_proc
    if f == 'fa':
        file_name = f'{ucsc_name}.fa.gz'
    elif f == 'gtf':
        file_name = f'genes/{ucsc_name}.ncbiRefSeq.gtf.gz'

    #download file
    cmd = f'wget https://hgdownload.soe.ucsc.edu/goldenPath/{ucsc_name}/bigZips/{file_name} -P {refD} -o {refD}/log.download_{f}.txt'
    call_return =0
    call_return = sub.call(cmd, shell=True)

    if call_return != 0:
        return call_return
    elif call_return == 0:
        file_name = file_name.split("/")[-1]
        sub.call(f'gzip -d {refD}/{file_name}', shell=True)
        return call_return



