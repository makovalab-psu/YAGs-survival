# python ../../assign_palindrome_rename.py CDY_CDYL_CDYL2.aligned_trimmed.fasta
# python ../../assign_palindrome_rename.py CDY_CDYL.aligned_trimmed.fasta
# python ../../assign_palindrome_rename.py CDY.aligned_trimmed.fasta


# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/CDY_CDYL_CDYL2.aligned_trimmed.renamed.fasta -T 12 -bb 1000 > run.log
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/CDY_CDYL.aligned_trimmed.renamed.fasta -T 12 -bb 1000 > run.log
# /Users/kxp5629/tools/iqtree-2.3.6-macOS/bin/iqtree2 -s renamed/CDY.aligned_trimmed.renamed.fasta -T 12 -bb 1000 > run.log


 Rscript ../../plot_tree.R renamed/CDY_CDYL_CDYL2.aligned_trimmed.renamed.fasta.treefile CDY_CDYL_CDYL2

 Rscript ../../plot_tree.R renamed/CDY_CDYL.aligned_trimmed.renamed.fasta.treefile CDY_CDYL

 Rscript ../../plot_tree.R renamed/CDY.aligned_trimmed.renamed.fasta.treefile CDY