# metatranscriptomics
Our GitHub Metatranscriptomics Repository offers tools for analyzing RNA transcripts in microbial communities, revealing their active roles in different environments.

Steps:
1-rRNA Removal: SortMeRNA is employed to filter out any residual rRNA sequences that were not eliminated during experimental procedures, thus refining the metatranscriptomic data to predominantly include non-rRNA transcripts. 

2-Quality Control: Raw sequencing data undergo quality assessment via FastQC and MultiQC. Cutadapt plays a crucial role in excising adaptor sequences from the raw data and quality is further augmented by Prinseq, which sieves out low-quality reads and sequences of low complexity. 

4-Host DNA removal: Subsequently, BMTagger is harnessed to remove any host-originating sequences, minimizing potential contamination. 

5-Taxonomy Assignmnet: The ensuing high-quality, rRNA-depleted reads are then channeled through the MetaPhlAn 4 pipeline to determine taxonomic classifications. 

6-Functional Annotation: Finally, to elucidate the functional attributes of the microbial community, HUMAnN 3 is applied, enabling the quantification of gene family and metabolic pathway abundacies.
