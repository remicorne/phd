from tqdm.contrib.concurrent import process_map

def parallel_process(executor, cases, description="Processing"):
    """Executes the operations performed by {executor} once per c   ase in cases in parralel

    Args:
        executor (function): functio
        cases (list(args)): different arguments to be used by executor

    Returns:
        [executor(case) for case in cases]
    """
    results = process_map(executor, cases, desc=description, chunksize=1)
    return results
    