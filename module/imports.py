import warnings
import pingouin as pg
import scipy
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import pandas as pd
from module.utils import *
from module.getters import *
from module.statistics import *
from module.figures import *
from module.histogram import *
from outliers import smirnov_grubbs as grubbs
import matplotlib.pyplot as plt
