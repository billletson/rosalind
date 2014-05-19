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
