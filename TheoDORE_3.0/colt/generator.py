"""Basic Logic for an INI based code generator"""
from abc import abstractmethod
from collections.abc import Mapping, MutableMapping
#
import re
#
from .slottedcls import slottedcls
from .configast import parse, IS_BLOCK, Entry


class GeneratorNavigator:
    """Contains all core logic to navigate
       through a generator output

       Sometimes its usefull to have access only to this logic
       without the whole `Generator` class
    """
    #
    seperator = "::"
    # is there a branching in the tree?
    is_branching_regex = re.compile(r"(?P<branch>[\d\w-]*)\((?P<node>.*)\)")
    # named tuple to store the pair
    Branch = slottedcls("Branch", ["branch", "node"])

    @classmethod
    def get_node_from_tree(cls, node_name, tree):
        """return node from a given tree"""
        return cls._get_node(node_name, tree)

    @classmethod
    def get_branching(cls, node):
        """check if node is a branching node"""
        branch = cls.is_branching_regex.match(node)
        if branch is None:
            return branch
        return cls.Branch(branch.group("branch"), branch.group("node"))

    @classmethod
    def join_keys(cls, parent, node):
        if parent is None or parent == "":
            return node
        if node is None or node == "":
            return parent
        return f"{parent}{cls.seperator}{node}"

    @staticmethod
    def join_case(branch, case):
        return f"{branch}({case})"

    @classmethod
    def rsplit_keys(cls, node):
        parent, sep, child = node.rpartition(cls.seperator)
        if sep != cls.seperator:
            parent = child
            child = None
        return parent, child

    @classmethod
    def _get_node(cls, node_name, tree):
        """Actual implementation of get_node """
        if node_name is None or node_name.strip() == "":
            return tree
        # loop over nodes to reach the particular node
        nodes = node_name.split(cls.seperator)
        for node in nodes:
            tree = cls._get_next_node(tree, node)
            if tree is None:
                return None
        return tree

    @classmethod
    def _get_next_node(cls, tree, node):
        """Get the next node of the current tree
           using the key to the node!

           Args:
                tree (object):
                    Python object that behaves like tree with
                    only dictionaries inside or objects that
                    behave like dictionaries

                node (string):
                    string to determine the next node
        """

        conditions = cls.get_branching(node)
        if conditions is None:
            return tree.get(node, None)

        node, case = conditions.branch, conditions.node
        tree = tree.get(node, None)
        if tree is not None:
            return tree.get(case, None)
        return tree


class Generator(GeneratorNavigator, Mapping):
    """Contains all core logic to generate
       code from a given config string/file

       Main use is the QuestionGenerator in colt,
       but it can also be used for other config based
       automatic code generators

       **IMPORTANT NOTE**:
       This works only if the objects generated by this method are `Mappings`
       unless they are `leaf nodes`, which can be anything.

       The tree generated by the basic generator can contain the following
       elements:

                 Node: basic container for the elements of our tree
                       basically is a `dict` and should behave like one
                       every node can contain several subnodes.

            Leaf Node: the actual elements one wants to store within
                       the tree, they are stored at the end

            Branching: a conditional branching. Similar to a standard `Node`,
                       but a branching contains a Leaf Node and then a dict of `Nodes`.

        A tree can be any of the elements above, but typically it is
        either a Branching or a Node
    """
    __slots__ = ('tree', '_keys')
    # please select leaf and branching type
    leafnode_type = type(None)
    branching_type = type(None)
    #
    node_type = dict
    #

    def __init__(self, treeconfig, *, comment=None):
        """Main Object to generate abstract tree from configstring

        Args:
            config(string):
                string that should be converted in the tree
        """
        if isinstance(treeconfig, Generator):
            self.tree = treeconfig.tree
            self._keys = treeconfig._keys
            return
        #
        if not isinstance(treeconfig, str):
            raise TypeError("Generator only accepts type string!")
        self.tree, self._keys = self._configstring_to_keys_and_tree(treeconfig, comment)

    @abstractmethod
    def leaf_from_string(self, entry, *, parent=None):
        """Create a leaf from an entry in the config file

        Args:
            entry, AstEntry
                contains main info for value

        Kwargs:
            parent (str):
                identifier of the parent node

        Returns:
            A leaf node

        Raises:
            ValueError:
                If the value cannot be parsed
        """

    def new_branching(self, name, *, leaf=None):
        """Create a new empty branching"""
        raise NotImplementedError("Branching not implemented in this tree, "
                                  "please implement the 'new_branching' method")

    @staticmethod
    def new_node(comment=None):
        """Create a new node of the tree"""
        return {}

    @staticmethod
    def tree_container(comment=None):
        """Create the container of the Tree"""
        return {}

    def get_node(self, node_name):
        """Parse down the abstraction tree to extract
           a particular node based on its block name
        """
        return self._get_node(node_name, self.tree)

    def add_element(self, name, line, *, comment=None, parentnode=None, overwrite=False):
        """add a single leaf node to the tree"""
        tree = self._get_subtree(parentnode)
        print(tree)
        if name in tree:
            if overwrite is False:
                raise KeyError(f"Node {name} already exists in {parentnode}")
        tree[name] = self.leaf_from_string(Entry(name, line, comment), parent=parentnode)

    def add_elements(self, configtree, *, parentnode=None, overwrite=True):
        """add elements to a particular node of the tree"""
        tree = self._get_subtree(parentnode)
        # check that treeconfig is correct type!
        subtree, keys = self._get_keys_and_subtree(configtree, parentnode=parentnode)
        self._update_keys(keys)
        # update subtree
        if overwrite is True:
            tree.update(subtree)
            return
        # overwrite it
        for key, item in subtree.items():
            if key not in tree:
                tree[key] = item

    def add_node(self, name, config, *, parentnode=None):
        """Add a new `node` with name `name` in given `parentnode`

        Args:
            name (str):
                name of the new node

            config (string, tree):
                the config specifying the given node

        Kwargs:
            parentnode (str):
                The name of the parent node, the new node should be created in

        Raises:
            ValueError:
                If a node already exists

        Example:
            >>> _question = "sampling = "
            >>> questions.add_node("software", {name: software.questions for name, software
                                                in cls._softwares.items()})
        """
        tree = self._get_subtree(parentnode)
        subtree, keys = self._get_keys_and_subtree(config, name=name, parentnode=parentnode)
        self._update_keys(keys)
        #
        if tree.get(name) is None:
            tree[name] = subtree
        else:
            raise ValueError(f"Node '{name}' in [{parentnode}] should not exist")

    def add_branching(self, leaf_name, branching_cases, *, parentnode=None):
        """Add a branching node inside `parentnode` with name `leaf_name`

        Args:
            leaf_name (str):
                name of leaf that should be replaced with a branching node

            branching_cases (dict):
                dictionary containing the cases that are implemented in the branching

        Kwargs:
            parentnode (str):
                name of the parent node

        Example:
            >>> questions.add_branching("sampling", {name: sampling.questions for name, sampling
                                                     in cls._sampling_methods.items()})
        """
        tree = self._get_subtree(parentnode)
        # generate the new branching, in case it exist return exisiting one
        tree[leaf_name] = self._new_branching_node(tree.get(leaf_name), leaf_name)
        # add branching to keys
        cases = self.join_keys(parentnode, leaf_name)
        self._update_keys(cases)
        # add cases!
        for case, config in branching_cases.items():
            # create name of real parent
            parent = self.join_keys(parentnode, self.join_case(leaf_name, case))
            subtree, keys = self._get_keys_and_subtree(config, parentnode=parent)
            self._update_keys(keys)
            tree[leaf_name][case] = subtree

    def __getitem__(self, key):
        """Mapping logic"""
        node = self.get_node(key)
        if node is None:
            raise KeyError(f"Node {key} does not exisit")
        return node

    def __len__(self):
        return len(self._keys)

    def __iter__(self):
        """return sorted keys just in case"""
        return iter(sorted(self._keys))

    def block_items(self):
        for key, value in self.items():
            if isinstance(value, self.node_type):
                yield key, value

    @staticmethod
    def _preprocess_string(string):
        """Basic Preprocessor"""
        return string

    def _configstring_to_keys_and_tree(self, string, comment):
        """transform a configstring to a tree object"""
        return self._generate_tree(string, comment)

    @classmethod
    def _get_parent_node(cls, node, tree):
        """Go iterative through the tree and get
           the final node, one over the selected one
           This is done in case the selected node does not
           exist and should be created.

           Args:
                node (str):
                    String to find the corresponding node in the tree

                tree (object):
                    python object that behaves like a dict of dicts
            Returns:
                final_node (str):
                    identifier of the final node
                tree (object):
                    tree object of the parent node
        """
        #
        nodes = node.split(cls.seperator)
        final_node = nodes[-1]
        #
        for nodename in nodes[:-1]:
            tree = cls._get_next_node(tree, nodename)
            if tree is None:
                return None, None
        return final_node, tree

    def _select_subnode(self, tree, nodeblock, *, comment=None):
        """Get a node from the tree creates:
           if it does not exist inside the parent node!

           iterative loop over the nodes till the selected one is reached
           in case the parent node is a branching, create a new one, if
           it is not done yet!

        """
        node, tree = self._get_parent_node(nodeblock, tree)
        if node is None:
            return None
        # if is not decission, create the new node as an dict
        conditions = self.get_branching(node)
        if conditions is None:
            if node in tree:
                block, _, _ = nodeblock.rpartition(self.seperator)
                raise KeyError(f"{block} already exists in {nodeblock}")
            tree[node] = self.new_node(comment=comment)
            return tree[node]
        #
        branch_name, node_name = conditions.branch, conditions.node
        branching = tree.get(branch_name, None)
        # insert new branching into tree
        branching = self._new_branching_node(branching, branch_name)
        if node_name not in branching:
            branching[node_name] = self.new_node(comment=comment)
        tree[branch_name] = branching
        return branching[node_name]

    def _new_branching_node(self, branching, branch_name):
        """
        Return a branching node, if the branching does not exist create it!

        Args:
            branching (obj):
                can be None, leaf node, branching

                    None: selected branching does not exist.
                          Create branching and a node
                          with name `node_name` inside the branching

               Leaf Node: selected branching does not exist, but a leaf node at that position
                          Create branching from that leaf node and init a node
                          with name `node_name` inside that branching

               Branching: selected branching does exist, and
                          an empty node will be created in that branching

        """
        if branching is None:
            branching = self.new_branching(branch_name)
        elif isinstance(branching, self.leafnode_type):
            branching = self.new_branching(branch_name, leaf=branching)
        elif isinstance(branching, self.branching_type):
            pass
        else:
            raise TypeError("Branching can only be typ: ",
                            f"'None', '{self.leafnode_type}', '{self.branching_type}'")
        #
        return branching

    @classmethod
    def _is_subblock(cls, block):
        """Check if block is further than one node away from
           the tree root """
        if any(key in block for key in (cls.seperator, '(', ')')):
            return True
        return False

    def _update_keys(self, keys):
        if isinstance(keys, set):
            self._keys |= keys
        elif isinstance(keys, str):
            self._keys.add(keys)
        else:
            TypeError(f"keys can only be set or str not {type(keys)}")

    def _get_keys(self, iterator, *, name=None, parentnode=None):
        keys = set()
        if name is not None:
            name = self.join_keys(parentnode, name)
            keys.add(name)
        else:
            name = parentnode
        keys |= set(self.join_keys(name, key) for key in iterator)
        return keys

    def _get_subtree(self, node_name):
        """get the node of a subtree at a given position

           important, it can not return branchings, only subtrees
        """
        subtree = self.get_node(node_name)
        # check that subtree is correct!
        if subtree is None:
            raise KeyError(f"Node {node_name} unknown")
        if not isinstance(subtree, self.node_type):
            raise ValueError(f"Node {node_name} has to be of type {self.node_type}")
        return subtree

    def _get_keys_and_subtree(self, configtree, *, name=None, parentnode=None):
        """update the tree keys and get a particular node"""
        if isinstance(configtree, Generator):
            keys = self._get_keys(configtree.keys(), name=name, parentnode=parentnode)
            subtree = configtree.tree
        elif isinstance(configtree, str):
            subtree, keys = self._configstring_to_keys_and_tree(configtree)
        else:
            raise TypeError("Generator only accepts type string or Generator!")
        return subtree, keys

    def _generate_tree(self, config, comment):
        """Generate a new tree from a configparser object

        Args:
            config (dict):
                linear dictionary of the corresponding config

        Returns:
            tree (object):
                parsed tree from the config file
        """
        # list to store keys in the created tree
        keys = set()
        # add starting value
        keys.add("")
        # create container for the tree
        maintree = self.tree_container(comment=comment)
        tree = maintree
        parent = None
        # parse defaults
        for entry in parse(config):
            if entry.value is IS_BLOCK:
                keys.add(entry.name)
                parent = entry.name
                tree = self._create_block(maintree, entry)
                continue
            tree[entry.name] = self.leaf_from_string(entry, parent=parent)
        return maintree, keys

    def _create_block(self, tree, entry):
        # get subsections
        return self._select_subnode(tree, entry.name, comment=entry.comment)


class BranchingNode(MutableMapping):
    """Basic Class that can be used as a Branching Node"""

    __slots__ = ('name', 'leaf', 'subnodes')

    def __init__(self, name, leaf, subnodes=None):
        self.name = name
        self.leaf = leaf
        if subnodes is None:
            subnodes = {}
        self.subnodes = subnodes

    def __getitem__(self, key):
        return self.subnodes[key]

    def __setitem__(self, key, value):
        self.subnodes[key] = value

    def __delitem__(self, key):
        del self.subnodes[key]

    def __iter__(self):
        return iter(self.subnodes)

    def __len__(self):
        return len(self.subnodes)

    def __str__(self):
        return (f"BranchingNode(name = {self.name},"
                f" leaf = {self.leaf}, subnodes = {self.subnodes}")

    def __repr__(self):
        return (f"BranchingNode(name = {self.name},"
                f" leaf = {self.leaf}, subnodes = {self.subnodes}")
