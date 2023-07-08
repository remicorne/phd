import warnings
import pingouin as pg
import scipy
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import pandas as pd



# The following functions are just here to be passed to the pd.corr() method, c and y are the two lists of values (columns) to be correlated
# This is a classic design pattern of which i forget the name again.
def isSignificant(stat_method_cb, pval_threshold=0.05):
    # As you can see here it will return a function. NOT CALL THE FUNCTION, but return it. The point here is to inject variables in the function that is returned.
    # when isSignificant(callback, pval_threshold) is called, it declare the anonymous function (lambda) passing the variables, and returns this declaration
    # this means that when the lambda function is called later on these will be 'harcoded' in the sense that they are no longer variables to be passed
    # Why? because if you look at the usage I pass it to pd.corr() to generate the mask in getPearsonCorrStats(). The thing is that pd.corr() calls the method it is given with (x, y) arguments
    # For the mask to be properly generated however, we also want to give a pval threshold to pass, as well as the statistical method that determines significance
    # But there is no way to pass this once you are in th pd.corr() context, so this is why you use this design pattern. It is meant to import variable into a context where they are not normally available
    return lambda x, y: stat_method_cb(x, y) >= pval_threshold


# The


def getPearson(x, y):
    return scipy.stats.pearsonr(x, y)


def getPearsonR(x, y):
    return getPearson(x, y).statistic


def getPearsonPValue(x, y):
    return getPearson(x, y).pvalue


@select_params
def getTukey(data, p_value_threshold):
    columns, *stats_data = pairwise_tukeyhsd(
        endog=data["value"], groups=data["treatment"], alpha=p_value_threshold
    )._results_table.data
    return pd.DataFrame(stats_data, columns=columns)


@select_params
def getOneWayAnova(data):
    F_value, p_value = scipy.stats.f_oneway(
        *[list(group_df["value"]) for treatment, group_df in data.groupby("treatment")]
    )
    # print(f'oneWAY_ANOVA F_value: {F_value}, p_value: {p_value}')
    return pd.DataFrame([[F_value, p_value]], columns=["F", "p_value"])


@select_params
def getTwoWayAnova(data, independant_vars):
    data[independant_vars] = data.apply(
        lambda x: [var in x["treatment"] for var in independant_vars],
        axis=1,
        result_type="expand",
    )
    return pg.anova(
        data=data,
        dv="value",
        between=independant_vars,
        detailed=True,
    ).round(3)


def doRawDfGrubbs(raw_df):
    result_list = []
    for group in raw_df.groupby("treatment"):
        result_list.append(
            grubbsTest(raw_df)
        )  # this func will loop through the df to feel grubbsTest
    return result_list


def grubbsTest(group_list):  # include the vairable type in name i.e. group_list series
    # from outliers import smirnov_grubbs as grubbs
    # x_ = grubbs.test(x, alpha=p_value)
    return



QUANTITATIVE_STAT_METHODS = {
    "twoway_anova": getTwoWayAnova,
    "oneway_anova": getOneWayAnova,
    "tukey": getTukey,
}


QUALITATIVE_STAT_METHODS = {"pearson": getPearson}


