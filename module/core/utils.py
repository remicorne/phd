from tqdm.contrib.concurrent import process_map
import itertools
from collections.abc import Iterable


def is_array_like(value):
    """Detect wheher a value is a non-string iterable

    Args:
        value (Any): Value to evaluate

    Returns:
        bool: True if value is a non-string iterable
    """
    return isinstance(value, Iterable) and not isinstance(value, str)


def call_case(case):
    """
    Calls the given case and returns its result.
    Used by parallel_process for complex calculation where passing case to a function is complex.

    Args:
        case (Callable): A callable object to be executed.

    Returns:
        Any: The result of the executed case.
    """
    return case()


def parallel_process(cases, executor=call_case, description="Processing"):
    """Executes the operations performed by {executor} once per case in cases in parralel

    Args:
        executor (function): function to be executed, must take a single argument and unpack if necessary
        cases (list(args)): different arguments to be used by executor. Case == [arg1, arg2, ...] for executor(case)

    Returns:
        List of results from executor
    """
    results = process_map(executor, cases, desc=description, chunksize=1)
    return results


def flatten(two_dimension_list):
    return list(itertools.chain.from_iterable(two_dimension_list))
