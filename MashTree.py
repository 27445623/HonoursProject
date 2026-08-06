from Bio import Phylo
tree = Phylo.read("mashtree.dnd", "newick")
print(tree)
Phylo.draw_ascii(tree)