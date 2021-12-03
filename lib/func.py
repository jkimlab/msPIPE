#! /usr/bin/env python3.6

import sys
import argparse
import configparser
import subprocess as sub

# ----------- system function ----------------

def printlog(logF, lines):
    O = open(logF, 'a')
    O.write(lines+'\n')
    print(lines)
    O.close()

def subproc(CMD, logF, log_proc):
    if logF == "stdout":
        printlog(log_proc, f'$ {CMD}')
        call = sub.call(CMD,  shell=True,universal_newlines=True)
    else:
        printlog(log_proc, f'$ {CMD}')
        call = sub.call(f'{CMD} 2> {logF}', shell=True, universal_newlines=True)

    if call == 0 :
        printlog(log_proc, '\t .. done\n')
        return 0
    else :
        printlog(log_proc, "\n*** There's a problem executing the upper command.")
        printlog(log_proc, f"*** check logfile > {logF}\n")
        return 1


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
    if f == 'fa':
        file_name = f'{ucsc_name}.fa.gz'
    elif f == 'gtf':
        file_name = f'genes/{ucsc_name}.ncbiRefSeq.gtf.gz'

    #download file
    cmd = f'wget https://hgdownload.soe.ucsc.edu/goldenPath/{ucsc_name}/bigZips/{file_name} -P {refD} -o {refD}/log.download_{f}.txt'
    call_return = sub.call(cmd, shell=True)

    if call_return != 0:
        return call_return
    elif call_return == 0:
        file_name = file_name.split("/")[-1]
        sub.call(f'gzip -d {refD}/{file_name}', shell=True)
        return call_return



