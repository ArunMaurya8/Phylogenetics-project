import numpy as np
import scipy.linalg as lin
import math


# Aim 1: store the tree.

class Tree:
    "Generic tree nodes."
    
    def __init__(self, name, children=None, parent=None, branch_length=None, sequence=None, nucleotide_matrix=None, vec_anc_lst = None): 
        self.name = name
        self.children = []
        if children != None:
            for child in children:
                self.add_children(child)
        self.parent = parent
        self.branch_length = branch_length
        self.sequence = sequence
        self.nucleotide_matrix = nucleotide_matrix
        self.vec_anc_lst = vec_anc_lst

    def add_children(self, child):
        self.children.append(child)

    def print_node(self):
        print('this is node ' + str(self.name) + ', sequence is ' + str(self.sequence))
        print('node ' + str(self.name) + ' has parent ' + str(self.parent) + ' and child(ren) ' + str(self.children) + ' with length ' + str(self.branch_length))
        print('vec_anc_lst:', self.vec_anc_lst)

# import data to the tree

with open('./data/table.dat') as t:
    table = t.readlines()
with open('./data/branchlength.dat') as l:
    lengths = l.readline()
with open('./data/msa.dat') as s:
    sequence = s.readlines()
length_lst = lengths.split(',')
pair_nodes = [[table[i].strip().split(',')[0], table[i].strip().split(',')[1]] for i in range(len(table))]
sequence_dir = {sequence[i].strip().split(' ')[0]:sequence[i].strip().split(' ')[1] for i in range(len(sequence))}
dict_base = {'A':[1,0,0,0], 'C':[0,1,0,0], 'G':[0,0,1,0], 'T':[0,0,0,1]}
node_dir = {}

for i in range(len(pair_nodes)):
    try:
        node_dir[pair_nodes[i][1]] = Tree(pair_nodes[i][1], parent=pair_nodes[i][0], branch_length=float(length_lst[i]), sequence=sequence_dir[pair_nodes[i][1]])
        node_dir[pair_nodes[i][1]].nucleotide_matrix = [dict_base[s] for s in node_dir[pair_nodes[i][1]].sequence]
    except KeyError:
        node_dir[pair_nodes[i][1]] = Tree(pair_nodes[i][1], parent=pair_nodes[i][0], branch_length=float(length_lst[i]))
for i in range(len(pair_nodes)):
    try:
        node_dir[pair_nodes[i][0]].add_children(pair_nodes[i][1])
    except KeyError:
        node_dir[pair_nodes[i][0]] = Tree(pair_nodes[i][0], children=[pair_nodes[i][1]])

# Aim 2: calculate the likelihood of the tree

# calculation for 1 site
def log_p_calculation(branch_lens, nucleotides):
    Q = np.array([[-0.5625,0.1875,0.1875,0.1875], [0.1875,-0.5625,0.1875,0.1875], [0.1875,0.1875,-0.5625,0.1875], [0.1875,0.1875,0.1875,-0.5625]])
    p_value_1 = lin.expm(Q * branch_lens[0])
    p_value_2 = lin.expm(Q * branch_lens[1])
    vec_anc_lst = []
    for s in range(len(nucleotides[0])):
        vec_val_1 = None
        vec_val_2 = None
        vec_anc = None
        vec_val_1 = np.matmul(p_value_1, nucleotides[0][s])
        vec_val_2 = np.matmul(p_value_2, nucleotides[1][s])
        vec_anc = [vec_val_1[i] * vec_val_2[i] for i in range(4)]
        vec_anc_lst.append(vec_anc)

    return vec_anc_lst
# calculate the vec_anc value for each ancestor including the internal ancestor
for key, value in node_dir.items():
    if len(value.children) != 0:
        if node_dir[value.children[0]].nucleotide_matrix != None:
            node_dir[key].vec_anc_lst = log_p_calculation([node_dir[value.children[0]].branch_length, node_dir[value.children[1]].branch_length], [node_dir[value.children[0]].nucleotide_matrix, node_dir[value.children[1]].nucleotide_matrix])
        elif node_dir[value.children[0]].nucleotide_matrix == None and node_dir[value.children[1]].nucleotide_matrix == None: 
            node_dir[key].vec_anc_lst = log_p_calculation([node_dir[value.children[0]].branch_length, node_dir[value.children[1]].branch_length], [node_dir[value.children[0]].vec_anc_lst, node_dir[value.children[1]].vec_anc_lst])
        elif node_dir[value.children[0]].nucleotide_matrix == None and node_dir[value.children[1]].nucleotide_matrix != None: 
            node_dir[key].vec_anc_lst = log_p_calculation([node_dir[value.children[0]].branch_length, node_dir[value.children[1]].branch_length], [node_dir[value.children[0]].vec_anc_lst, node_dir[value.children[1]].nucleotide_matrix])

# calculate the log likelihood for the root node            
vec_root = node_dir['6'].vec_anc_lst
log_tree = 0
for i in range(len(vec_root)):
    log_tree += np.log(np.matmul(vec_root[i], [0.25,0.25,0.25,0.25]))  
print('The log likelihood of the tree is: ', log_tree)      