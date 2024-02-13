#!/usr/bin/env python

import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats
from Bio import SeqIO
from statistics import mean
from operator import itemgetter


test_dict = {'gfg': 1, 'is': 4, 'best': 6, 'for': 7, 'geeks': 3}

res = dict(sorted(test_dict.items(), key=itemgetter(1), reverse=True)[:int(3)])

print(res)