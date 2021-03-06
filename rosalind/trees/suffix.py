"""
Implementation of Ukkonen's algorithm, a linear time Suffix Tree construction algorithm
See: http://www.stanford.edu/~mjkay/gusfield.pdf
"""


class SuffixTree:
    def __init__(self, word, end_of_word="$"):
        if word[-1] != end_of_word:
            self.word = word + end_of_word
        else:
            self.word = word
        self._e = SuffixEndPoint()
        self.nodes = [SuffixNode(0, 0, downstream=[1]), SuffixNode(0, self._e, upstream=0)]
        rule3 = 0
        for i in xrange(1, len(self.word)):
            inserted = None
            active_node = 0
            for j in xrange(rule3, i + 1):
                active_node, rule3, inserted = self._extend(j, i, active_node, inserted)
                if rule3:
                    break
            self._e.increment()

    def __repr__(self):
        return "\n".join(["String: " + self.word] + [repr(x) for x in self.nodes])

    def _extend(self, j, i, active_node=0, inserted=None):
        """
        Extends and returns:
        The node to which the suffix was added (if added inside an edge, the node above the edge)
        If rule three was invoked, meaning we can skip the rest of the step
        If a new node was created from splitting an edge, what that node is, allowing for creation of a suffix link
        """

        # Node 1 is always the full string, can increment its length each time
        if j == 0:
            return 1, 0, None

        # Walk back up until hit the root or a suffix link. Record how far we traveled.
        active_node, back_distance = self._backtrack(active_node)

        # Skip count trick, recording j first
        j0 = j
        if back_distance >= 0:
            j = i - back_distance
        end_node, end_index = self._skip_count(j, i, active_node)

        #Figure out how to extend
        if not self.nodes[end_node].downstream and end_index + 1 == self.nodes[end_node].edge_length:
            # Rule 1 - ending at a leaf (not along the edge leading to a leaf)
            return end_node, 0, None
        elif end_index + 1 == self.nodes[end_node].edge_length:
            for x in self.nodes[end_node].downstream:
                if self.word[self.nodes[x].edge_start] == self.word[i]:
                    # Rule 3 - Letter to insert is the first letter of an edge leading to a child node
                    return x, j0 - 1, None
            else:
                # Rule 2, with branch at current node (current node has children, but none start with letter to
                # be inserted)
                self.nodes[end_node].downstream.append(len(self.nodes))
                self.nodes.append(SuffixNode(i, self._e, end_node))
                return len(self.nodes) - 1, 0, None
        else:
            if self.word[self.nodes[end_node].edge_start + end_index + 1] != self.word[i]:
                # Rule 2, current ends within an edge, next letter on the edge is not the letter to be inserted
                above = self.nodes[end_node].upstream
                mid = len(self.nodes)
                self.nodes[above].downstream.remove(end_node)
                self.nodes[above].downstream.append(mid)
                self.nodes[end_node].upstream = mid
                self.nodes.append(SuffixNode(self.nodes[end_node].edge_start,
                                             self.nodes[end_node].edge_start + end_index + 1, above,
                                             [end_node, mid + 1]))
                self.nodes.append(SuffixNode(i, self._e, mid))
                self.nodes[end_node].edge_start += end_index + 1
                if inserted is not None:
                    self.nodes[inserted].suffix_link = mid
                return mid, False, mid
            else:
                # Rule 3, letter to be inserted is the next along the current edge
                return end_node, j0 - 1, None

    def _backtrack(self, active_node):
        """
        Travel back up the tree until we hit the root (stop)  or a suffix_link (traverse). Return the new active_node,
        and how many characters we backtracked (return -1 if hit the root).
        Arguments: int active_node
        Returns: int, int
        """
        if active_node == 0:
            return 0, -1
        back_distance = 0
        while True:
            if self.nodes[active_node].suffix_link is not None:
                return self.nodes[active_node].suffix_link, back_distance
            elif self.nodes[active_node].upstream == 0:
                return 0, -1
            else:
                back_distance += self.nodes[active_node].edge_length
                active_node = self.nodes[active_node].upstream

    def _skip_count(self, j, i, active_node):
        """
        Traverse the tree downwards, only checking the first letter of an edge since we are guaranteed to have the
        string we are searching for in the tree (as such, we take indicies into the stored word rather than a string).
        Returns the new active_node, and the index of where we are along its edge.
        Arguments: int j, int i, int active_node
        Returns int, int
        """
        distance = i - j
        while True:
            for node in self.nodes[active_node].downstream:
                if self.word[self.nodes[node].edge_start] == self.word[j]:
                    if self.nodes[node].edge_length < distance:
                        distance -= self.nodes[node].edge_length
                        j += self.nodes[node].edge_length
                        active_node = node
                        break
                    elif self.nodes[node].edge_length == distance:
                        return node, self.nodes[node].edge_length - 1
                    else:
                        return node, distance - 1
            else:
                return active_node, self.nodes[active_node].edge_length - 1

    def find_path(self, string, active_node=0):
        """
        Traverse the tree downwards, checking each letter in turn. Return the complete path, the ending node, and the
        index of where we are along its edge. If the string is not in the tree, returns (None, None, None).
        Arguments: str string, int active_node
        Returns int[], int, int
        """
        if not string:
            return [0], 0, -1
        edge_index = self.nodes[active_node].edge_length
        path = [active_node]
        for letter in string:
            if edge_index >= self.nodes[active_node].edge_length:
                for x in self.nodes[active_node].downstream:
                    if letter == self.word[self.nodes[x].edge_start]:
                        active_node = x
                        edge_index = 1
                        path.append(x)
                        break
                else:
                    return None, None, None
            else:
                if letter == self.word[self.nodes[active_node].edge_start + edge_index]:
                    edge_index += 1
                else:
                    return None, None, None
        return path, active_node, edge_index - 1

    def _count_leaves(self, node):
        """
        Counts the number of leaves below a node, equivalent to the occurences of the substring leading to that node.
        Arguments: int
        Returns: int
        """
        if self.nodes[node].downstream:
            return sum([self._count_leaves(x) for x in self.nodes[node].downstream])
        else:
            return 1

    def substring_count(self, substring):
        """
        Counts the occurences of a substring in the tree
        Arguments: str
        Returns: int
        """
        path, node, index = self.find_path(substring)
        return self._count_leaves(node)

    def substring_at_node(self, node):
        """
        Reads the substring that leads from the root to the specified node.
        Arguments: int
        Returns: str
        """
        substring = ""
        while node != 0:
            # subtraction here is a hack because end can be an int or object
            substring = self.word[self.nodes[node].edge_start:self.nodes[node].edge_end - 0] + substring
            node = self.nodes[node].upstream
        return substring

    def longest_occuring_k_times(self, k):
        """
        Finds the longest substring that occurs at least k times in the string. If a tie, returns any.
        Arguments: int
        Returns: str
        """
        potentials = [self.substring_at_node(node) for node in xrange(len(self.nodes)) if self._count_leaves(node) >= k]
        return max(potentials, key=len)

    def list_of_edges(self):
        """
        Returns a list of all the edges in the tree. These are not full substrings, just the edge labels.
        Arguments:
        Returns: [str]
        """
        return [self.word[node.edge_start:node.edge_end - 0] for node in self.nodes[1:]]

    def list_of_substrings(self):
        """
        Returns a list of all the substrings in the tree.
        Arguments:
        Returns: [str]
        """
        return [self.substring_at_node(node) for node in xrange(1, len(self.nodes))]


class SuffixNode:
    def __init__(self, edge_start, edge_end, upstream=None, downstream=[], suffix_link=None):
        self.edge_start = edge_start
        self.edge_end = edge_end
        self.upstream = upstream
        self.downstream = downstream
        self.suffix_link = suffix_link

    @property
    def edge_length(self):
        return self.edge_end - self.edge_start

    def __repr__(self):
        return "Start:%s, End:%s, Upstream:%s, Downstream:%s, Suffix Link:%s" % \
               (self.edge_start, self.edge_end, self.upstream, self.downstream, self.suffix_link)


class SuffixEndPoint:
    def __init__(self, e=1):
        self.e = e

    def __sub__(self, other):
        return self.e - other

    def increment(self):
        self.e += 1

    def __repr__(self):
        return "e" + str(self.e)
