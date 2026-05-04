# python  ../assign_palindrome_rename.py DAZ_DAZL_RBD.aligmnent.20240924.fasta
# python  ../assign_palindrome_rename.py DAZ_DAZL_repeats.alignment.20240924.fasta

# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/DAZ_DAZL_RBD.aligmnent.20240924.renamed.fasta -T 12 -bb 1000 > run.log
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/DAZ_DAZL_repeats.alignment.20240924.renamed.fasta -T 12 -bb 1000 >> run.log

Rscript ../plot_tree.R renamed/DAZ_DAZL_RBD.aligmnent.20240924.renamed.fasta.treefile DAZ
# Rscript ../plot_tree.R renamed/DAZ_DAZL_repeats.alignment.20240924.renamed.fasta.treefile DAZ SymSyn_chr10 40