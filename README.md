# Phylogenetics project
Aim: Build a python module to calculate the likelihood of a phylogenetic tree (guided by Prof Nicolas Salamin)

The script is to store the information about the tree and calculate the log-likelihood of the tree from the given information. The data consists of following three files:

1. an MSA in fasta (msa.dat)
2. a tree in table format (table.dat)
3. a list of branch lengths (branchlength.dat)

The task is to:

Step1: 
- create a tree object to store the phylogeny
- import the parent-child table and the branch lengths vector in python
- populate the tree object with each node in the table

(A tree is a list of nodes. Nodes have attributes like branch lengths, parent, children etc.)
![image](https://github.com/user-attachments/assets/09409a9d-3a45-45f1-a8f5-864d6810372b)


Step2: Use dynamic programming to calculate the probability of the alignment on the tree.

![image](https://github.com/user-attachments/assets/0b5993e6-b651-4662-b5c3-3c13f7239f5e)


To calculate the log-likelihood of the known tree mentioned, run the command:

```python3 script.py```

