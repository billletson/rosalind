def n_connected_subgraphs(nodes, edges):
    """
    Given a set of nodes and undirected edges, assigns each node
    to a subgraph and returns the number of connected subgraphs
    which are not connected to each other
    Arguments: int nodes, (int,int)[] edges
    Returns: int
    """
    assignments = {}
    for n in nodes:
        assignments[n] = n
    for edge in edges:
        assign_to = assignments[edge[0]]
        assign_from = assignments[edge[1]]
        for a in assignments.items():
            if a[1] == assign_from:
                assignments[a[0]] = assign_to
    return len(set(assignments.values()))


def unrooted_internal_from_leaves(leaves):
    return leaves - 2


def unrooted_leaves_from_internal(internal):
    return internal + 2


def overlap_graph(dnas, o):
    """
    Compute overlap graph given a list of Sequence objects and a
    prefix/suffix length o. Returns list of tuples of directed connections
    Arguments: Sequence[], int
    Returns: [(str,str)]
    """
    return [(x.name, y.name) for x in dnas
                            for y in dnas
                                 if (x.name != y.name or x.sequence != y.sequence)
                                     and x.sequence[-o:] == y.sequence[:o]]