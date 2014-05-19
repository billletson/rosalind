class NewickTree:
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
    def _separate_weight(cls, text):
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
    def _find_common_ancestor_path(cls, first_path, second_path):
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

