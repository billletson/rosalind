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


class Trie:
    def __init__(self, words):
        self.max_node = 0
        self.children = [[]]
        self.edge_labels = {}
        for word in words:
            self.insert(word)

    def insert(self, word):
        current_node = 0
        for letter in word:
            for i in self.children[current_node]:
                if letter == self.edge_labels[(current_node, i)]:
                    current_node = i
                    break
            else:
                self.max_node += 1
                self.children[current_node].append(self.max_node)
                self.children.append([])
                self.edge_labels[(current_node, self.max_node)] = letter
                current_node = self.max_node

    def edges_as_strings(self, one_index=False):
        if one_index:
            offset = 1
        else:
            offset = 0
        strings = [" ".join([str(x[0][0] + offset), str(x[0][1] + offset), x[1]]) for x in self.edge_labels.items()]
        return strings


class Tree:
    def __init__(self, newick):
        self.root = self._parse_to_node(newick)

    def _parse_to_node(self, newick):
        if not newick:
            return Node(newick)
        if newick[-1] == ";":
            newick = newick[:-1]
        if newick[0] != "(":
            name, weight = self._separate_weight(newick)
            return Node(name, weight)
        for i in xrange(len(newick)-1, 0, -1):
            if newick[i] == ")":
                name, weight = self._separate_weight(newick[i + 1:])
                node = Node(name, weight)
                remainder = newick[1:i]
                children = []
                paren_count = 0
                last_split = -1
                for j in xrange(len(remainder)):
                    if remainder[j] == "(":
                        paren_count += 1
                    elif remainder[j] == ")":
                        paren_count -= 1
                    elif remainder[j] == "," and paren_count == 0:
                        children.append(remainder[last_split + 1:j])
                        last_split = j
                children.append(remainder[last_split + 1:])
                node.children = [self._parse_to_node(child) for child in children]
                return node

    @classmethod
    def _separate_weight(self, text):
        if ":" in text:
            name, weight = text.split(":")
            weight = float(weight)
        else:
            name = text
            weight = 1
        return name, weight

    def find_distance(self, first, second):
        first_path = self._find_path(first)
        second_path = self._find_path(second)
        common_path = self._find_common_ancestor_path(first_path, second_path)
        return self._distance_from_root_by_path(first_path) + self._distance_from_root_by_path(second_path) \
               - 2 * self._distance_from_root_by_path(common_path)

    def _distance_from_root_by_path(self, path):
        distance = 0
        current_node = self.root
        for branch in path:
            current_node = current_node.children[branch]
            distance += current_node.distance
        return distance

    @classmethod
    def _find_common_ancestor_path(self, first_path, second_path):
        common_path = []
        for i in xrange(min(len(first_path), len(second_path))):
            if first_path[i] == second_path[i]:
                common_path.append(first_path[i])
            else:
                break
        return common_path

    def nontrivial_characters(self):
        node_names = sorted(self.root.node_names())
        nodes = self.root.nodes_with_children()
        nodes.remove(self.root)
        name_characters = [node.node_names() for node in nodes]
        binary_characters = ["".join(["1" if name in nc else "0" for name in node_names]) for nc in name_characters]
        return binary_characters

    def _find_path(self, name):
        return self.root.find_path(name)


class Node:
    def __init__(self, name=None, distance=1):
        self.name = name
        self.distance = distance
        self.children = []

    def find_path(self, name):
        if self.name == name:
            return []
        for i in xrange(len(self.children)):
            path = self.children[i].find_path(name)
            if path is not None:
                return [i] + path
        else:
            return None

    def nodes_with_children(self):
        ret = []
        if self.children:
            ret.append(self)
            for child in self.children:
                ret += child.nodes_with_children()
        return ret

    def node_names(self):
        ret = []
        if self.name != "":
            ret.append(self.name)
        if self.children:
            for child in self.children:
                ret += child.node_names()
        return ret


def debruijn(dnas):
    adjacencies = []
    for dna in dnas:
        adjacencies.append((dna.sequence[:-1], dna.sequence[1:]))
        rc = dna.reverse_complement()
        adjacencies.append((rc.sequence[:-1], rc.sequence[1:]))
    return set(adjacencies)















