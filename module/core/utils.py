import sys
from tqdm.contrib.concurrent import process_map
from itertools import chain
from functools import partial
from collections.abc import Iterable
from tqdm import tqdm
import math


def is_array_like(value):
    """Detect wheher a value is a non-string iterable

    Args:
        value (Any): Value to evaluate

    Returns:
        bool: True if value is a non-string iterable
    """
    return (
        isinstance(value, Iterable)
        and not isinstance(value, str)
        and not isinstance(value, dict)
    )


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


def execute_cases(cases, executor=call_case):
    """Used to optimise paralel and reduce ovehead of creating one process per case"""
    return [executor(case) for case in cases]


def parallel_process(
    cases, executor=call_case, description="Processing", optimize=False
):
    """Executes the operations performed by {executor} once per case in cases in parralel.
    With optimization, it is assumed each case requires minimal computation and that the overhead of creating processes is significant.
    The cases are thus split into 20 batches and processed in parallel, with a fallback to consecutive execution if there are less than 20 cases.

    Args:
        executor (function): function to be executed, must take a single argument and unpack if necessary
        cases (list(args)): different arguments to be used by executor. Case == [arg1, arg2, ...] for executor(case)
        optimize (bool): If True, the cases are split into 20 batches and processed in parallel.
    Returns:
        List of results from executor
    """

    run_synchronous = sys.gettrace() is not None  # Debug mode
    is_batched = False

    if optimize:
        if len(cases) > 20:
            batch_size = math.ceil(
                len(cases) / 20
            )  # Probably possible to have 20 processes
            cases = [
                cases[i : i + batch_size] for i in range(0, len(cases), batch_size)
            ]
            executor = partial(execute_cases, executor=executor)
            description = description + f" (batches of {batch_size} elements)"
            is_batched = True
        else:
            run_synchronous = True

    results = (
        [executor(batch) for batch in tqdm(cases, desc=description)]
        if run_synchronous  # Remove paralel processing for debugging
        else process_map(
            executor,
            cases,
            desc=description,
            chunksize=1,
        )
    )
    return list(chain(*results)) if is_batched else results


def strtobool(val: str) -> int:
    """
    Convert a string to a boolean represented as 1 (true) or 0 (false).

    Accepted true values are: 'y', 'yes', 't', 'true', 'on', '1'
    Accepted false values are: 'n', 'no', 'f', 'false', 'off', '0'

    Raises ValueError if 'val' is anything else.
    """
    val = val.strip().lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return 1
    if val in ("n", "no", "f", "false", "off", "0"):
        return 0
    raise ValueError(f"invalid truth value {val!r}")
