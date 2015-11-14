"""
Function: Try to find folder with shot data. If not listed,
it needs to be written at line 22.
Author:
Cassio Amador (cassioamador at yahoo.com.br)
"""

from os import path, pardir, makedirs


def set_folder(tipo):
    """Sets data folder to 'parent directory + SHOTS + data type'"""
    loc = path.join(pardir, 'SHOTS')
    makedirs(loc,exist_ok=True)
    #'data_type' is the subfolder, depending on data type. default is a shot.
    data_type = 'shots'
    if tipo == 'mirror':
        data_type = 'tests'
    elif 'test' in tipo:
        data_type = 'tests'
    elif tipo == 'clean':
        data_type = 'cleaning_plasma'
    folder = path.join(loc, data_type)
    makedirs(folder, exist_ok = True)
    return path.join(loc, data_type)

def check_prof_folder(shot):
    loc = path.join(pardir, "PROC_FILES")
    prof_folder = path.join(loc, "%s" % shot)
    makedirs(loc, exist_ok = True)
    makedirs(prof_folder, exist_ok = True)
    return prof_folder
