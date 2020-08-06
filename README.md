# Bremia-lactucae-WY-effector-characterization

- This repository contains scripts to analyze oomycete genomes for WY effector prediction and characterization, specifically focusing on the causal agent of lettuce downy mildew, *Bremia lactucae*.

### Paper on BioRxiv: [Effector prediction and characterization in the oomycete pathogen Bremia lactucae reveal host-recognized WY domain proteins that lack the canonical RXLR motif](https://www.biorxiv.org/content/10.1101/679787v2.supplementary-material)

*Kelsey Wood, Munir Nur, Juliana Gil, Kyle Fletcher, Kim Lakeman, Ayumi Gothberg, Tina Khuu, Jennifer Kopetzky, Archana Pandya, Mathieu Pel, Richard Michelmore*


### [Supplemental motif table with all secreted genes, annotated with GeneID, sequence, species, and RXLR/EER/WY motif presence](https://github.com/mjnur/Bremia-lactucae-WY-effector-characterization/blob/master/motif_counting/20200805_Supplemental_motif_category_table.csv)

Notes on the pipeline: 
  - RXLR string searches were done on non-redundant secretomes with cleaved sequences, yes the CSV above contains the full sequences for completeness. 
  - The [ORFs](https://github.com/mjnur/Bremia-lactucae-WY-effector-characterization/tree/master/ORFs) folder contains the starting sequences of our pipeline, which were used as input to SignalP v4.1 (sensitive mode) to acquire secretomes and cleaved sequences and CD-HIT to reduce sequence redundancy
