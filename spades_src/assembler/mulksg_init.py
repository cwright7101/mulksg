#!/usr/bin/env python
 
############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
from os.path import abspath, dirname, realpath, join, isfile

source_dirs = ["", "truspades", "common"]

# developers configuration
mulksg_home = abspath(dirname(realpath(__file__)))
bin_home = join(mulksg_home, 'bin')
python_modules_home = join(mulksg_home, 'src')
ext_python_modules_home = join(mulksg_home, 'ext', 'src', 'python_libs')
mulksg_version = ''


def init():
    global mulksg_home
    global bin_home
    global python_modules_home
    global mulksg_version
    global ext_python_modules_home

    # users configuration (mulksg_init.py and mulksg binary are in the same directory)
    if isfile(os.path.join(mulksg_home, 'mulksg-core')):
        install_prefix = dirname(mulksg_home)
        bin_home = join(install_prefix, 'bin')
        mulksg_home = join(install_prefix, 'share', 'spades')
        python_modules_home = mulksg_home
        ext_python_modules_home = mulksg_home

    for dir in source_dirs:
        sys.path.append(join(python_modules_home, 'spades_pipeline', dir))

    mulksg_version = open(join(mulksg_home, 'VERSION'), 'r').readline().strip()


if __name__ == '__main__':
    mulksg_py_path = join(dirname(realpath(__file__)), 'mulksg.py')
    sys.stderr.write('Please use ' + mulksg_py_path + ' for running MULKSG genome assembler\n')
