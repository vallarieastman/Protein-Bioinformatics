#!/usr/bin/env md-311
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 25 12:26:15 2025

@author: val
"""

import os
import sys
import pandas as pd
import numpy as np
import requests
from pprint import pprint

from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split, KFold, \
    GroupKFold
from scipy.stats import pearsonr, mode
import matplotlib.pyplot as plt
import seaborn as sns

# Import data from protein sequence embedding repository.
sys.path.append("/home/stephan/protein-sequence-embedding-iclr2019") #replace this path with YOUR path based on where you downloaded the github repo
from src.alphabets import Uniprot21
from torch.nn.utils.rnn import PackedSequence
from src.utils import pack_sequences, unpack_sequences

import torch
%matplotlib inline