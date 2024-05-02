from tqdm.contrib.concurrent import process_map

def parallel_process(executor, cases, description="Processing"):
    """Executes the operations performed by {executor} once per case in cases in parralel

    Args:
        executor (function): function to be executed, must take a single argument and unpack if necessary
        cases (list(args)): different arguments to be used by executor. Case == [arg1, arg2, ...] for executor(case)

    Returns:
        List of results from executor
    """
    results = process_map(executor, cases, desc=description, chunksize=1)
    return results
    