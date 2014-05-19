def find_reversals(original, target):
    """
    Find the minimum number of reversals required to transform the original into the target, and return the reversals.
    See paper at http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.91.3123&rep=rep1&type=pdf
    Arguments: str/list original, str/list target
    Returns: (int, int)[]
    """
    if original == target:
        return []
    order = [target.index(x) for x in original]
    return list(reversed(_recurse_find_reversals(order)))


def _recurse_find_reversals(order):
    """
    Find the minimum number of reversals required to transform the original into the target, and return the reversals.
    This function works recursively, testing each possible reversal and seeing the minimum number of reversals needed
    after completing that reversal. Returns the reversals as tuples of the two points to swap at.
    Arguments: int[]
    Returns: (int, int)[]
    """
    breakpoints, negative_strip = _find_breakpoints(order)
    if not breakpoints:
        return []
    reversals = []
    for bp in breakpoints:
        if bp == -1:
            swaps = [order.index(0)]
        elif order[bp] == len(order) - 1:
            swaps = [order.index(order[bp] - 1)]
        elif order[bp] == 0:
            swaps = [order.index(order[bp] + 1)]
        else:
            swaps = [order.index(order[bp] - 1), order.index(order[bp] + 1)]
        for swap in swaps:
            candidate_order = order[:bp+1] + list(reversed(order[bp+1:swap+1])) + order[swap+1:]
            candidate_breakpoints, negative_strip = _find_breakpoints(candidate_order)
            reduction = len(breakpoints) - len(candidate_breakpoints)
            if len(candidate_breakpoints) == 0:
                return [(bp + 1, swap)]
            if reduction > 1 or (reduction == 1 and negative_strip):
                reversals.append(_recurse_find_reversals(candidate_order) + [(bp + 1, swap)])
    if reversals:
        return min(reversals, key=len)
    else:
        return [0] * (len(order) + 1)


def _find_breakpoints(order):
    """
    Find the breakpoints in an series of indicies. A breakpoint occurs where two adjacent positions are not
    consecutively increasing or decreasing. If two or more adjacent positions are consecutively decreasing, this is
    called a negative strip. If this exists, return True, otherwise return False.
    Arguments: int[] order
    Returns: int[], bool
    """
    breakpoints = []
    negative_strip = False
    if order[0] != 0:
        breakpoints.append(-1)
    for i in xrange(len(order) - 1):
        if abs(order[i] - order[i + 1]) != 1:
            breakpoints.append(i)
        elif order[i] - order[i + 1] == 1:
            negative_strip = True
    return breakpoints, negative_strip

