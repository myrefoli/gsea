# gsea
A Python script that performs genome set enrichment analysis (GSEA) given a preprocessed expression file, sample file, and KEGG file.

To give some biological context for what GSEA is and why scientists use GSEA, we'll start at the cell level. Each cell in the human body has the same genetic code, yet performs a diverse array of different specialized functions (i.e. a liver cell is not the same as a lung cell). 

This is possible because of the central dogma of biology: we know that the genetic code (DNA) is used to transcribe (or “express”) genes to RNA transcripts (messenger RNA or mRNA) that carry information about which proteins are made in a cell. 

Most mRNA then goes on to the ribosome to be translated into proteins. But given the availability of the mRNA as intermediate step, gene expression assays can be used to understand cell and/or tissue function, by assessing which genes are being transcribed in the cell at a given time. Gene expression can also provide insight into which genes may be involved in disease processes.

The data for this specific project was drawn from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25628.
