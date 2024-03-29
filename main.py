######## ONLY IMPORT THE USER ACCESSIBLE FUNCTIONS HERE ##########
 # REMI this need sto have outliers and stats same methong as quant does
from module.quantitative import (
    quantitativeHistogram,
    percentageVehiclesFig,
    quantitativeSummary,
    justStats
)
from module.correlogram import correlogram, corrSelector, plot_correlation
from module.utils import initiateFileSystem
from module.metadata import saveMetadata
from module.utils import subselectDf
from module.getters import (
    getCompoundAndRatiosDf,
    getAggregateStatsDf,
    getRegionSubclassification,
    getQuantitativeStats,
)
from module.graph_theory import network, networkDegreeDistribution

#
from module.histogram import  buildHistogram
from module.pca import pca

######## INIT ##########
# Start by checking filesystem has all the folders necessary for read/write operations (cache) or create them otherwise
initiateFileSystem()
