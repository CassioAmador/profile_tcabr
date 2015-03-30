"""
Function: Try to find folder with shot data. If not listed,
it needs to be written at line 22.
Author:
Cassio Amador (cassioamador at yahoo.com.br)
TODO: read info from '.ini' file, in home dir.
"""

from os import getcwd, path


def find_folder(tipo):
    """Automatically finds the folder with shot data. If not listed,
    it needs to be written at line 22, in shots_folder.py"""
    #'loc' is the folder where the data is/will be stored.
    cur_dir = getcwd()
    if 'cassio' in cur_dir:
        if "Users" in cur_dir:
            loc = '/Users/cassio/fisica/TCABR/dados'
        else:
            loc = '/home/cassio/fisica/Reflectometria/TCABR/dados'
    elif 'GRS' in cur_dir or 'IPFNIST' in cur_dir:
        loc = '/home/GRS/TCABR/data'
    else:
        # User defined
        loc = path.join(getcwd(), '..')
    #'data_type' is the subfolder, depending on data type.
    data_type = 'shots'
    if tipo == 'mirror':
        data_type = 'tests'
    elif 'test' in tipo:
        data_type = 'tests'
    elif tipo == 'clean':
        data_type = 'cleaning_plasma'
    return path.join(loc, data_type)
