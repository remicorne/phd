######## ONLY IMPORT THE USER ACCESSIBLE FUNCTIONS HERE ##########
from module.headtwitch import (
    headTwitchHistogram,
)  # REMI this need sto have outliers and stats same methong as quant does
from module.quantitative import (
    quantitativeHistogram,
    percentageVehiclesFig,
    quantitativeSummary,
)
from module.correlogram import correlogram
from module.utils import initiateFileSystem
from module.metadata import saveMetadata
from module.utils import subselectDf
from module.getters import (
    getCompoundAndRatiosDf,
    getAggregateStatsDf,
    getHeadTwitchDf,
    getRegionSubclassification,
    getQuantitativeStats,
)

#
from module.histogram import buildHeadTwitchHistogramData, buildHistogram


######## INIT ##########
# Start by checking filesystem has all the folders necessary for read/write operations (cache) or create them otherwise
initiateFileSystem()
