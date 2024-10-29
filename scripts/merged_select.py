import subprocess
import argparse
import os, sys, re
from pathlib import Path
import gzip, statistics
import fnmatch
import numpy as np
import scipy.stats as st
import pandas as pd
import math
from itertools import islice

fin = str(sys.argv[1])
fout = str(sys.argv[2])

with open(fout, 'w') as fileout:
    with open(fin, 'r') as filein:
        while True:
            lines_gen = list(islice(filein, 4))
            if not lines_gen:
                break
            # process lines_gen
            rlength = int(lines_gen[0].rstrip().split()[1].split('_')[1]) + int(lines_gen[0].rstrip().split()[1].split('_')[2])
            if rlength >= 150:
                for line in lines_gen:
                    if line.rstrip()[0] != '@' and line.rstrip()[0] != '+':
                        fileout.write(line.rstrip()[10:-10]+'\n')
                    else:
                        fileout.write(line.rstrip()+'\n')
                
