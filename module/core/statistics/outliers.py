import numpy as np
from outliers import smirnov_grubbs


def grubbs_test(values, p_value_threshold):
    """
    Takes a list of values on which to perform the test and returns outliers
    """
    if len(values) < 3:
        return []
    return smirnov_grubbs.two_sided_test_outliers(values, alpha=p_value_threshold)


def iqr_test(values, k):
    if len(values) < 3:
        return []
    q1 = np.percentile(values, 25)
    q3 = np.percentile(values, 75)
    iqr = q3 - q1
    lower_bound = q1 - k * iqr
    upper_bound = q3 + k * iqr
    outliers = [x for x in values if x < lower_bound or x > upper_bound]
    outliers = sorted(
        outliers, key=lambda x: abs(x - (q1 if x < lower_bound else q3)), reverse=True
    )
    return outliers


def get_outliers(values, method, p_value_threshold):
    if method not in OUTLIER_TESTS:
        raise ValueError(f"Unknown outlier method: {method}")
    return OUTLIER_TESTS[method](values, p_value_threshold)


OUTLIER_TESTS = {"grubbs": grubbs_test, "iqr": iqr_test}
