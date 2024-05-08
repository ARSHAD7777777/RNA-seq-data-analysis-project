# Introduction
In this project, we will look at the changes in gene expression in Phaeodactylum tricornutum depending on the nutrient availability in the environment.
This is done by using data from RNA-seq. The data is pre-processed and analyzed in R using EdgeR, Bioconductor, and Limma to draw insights into various genes' fold changes. 
P. tricornutum is a diatom and often serves as a model system for marine microalgae.
One of the macronutrients needed to sustain the growth of marine algae populations is phosphor-
ous (P). In biological molecules, P is often represented in the form of covalent-bound phosphate
(PO3−4 ). However, inorganic phosphate (Pi) is also present in the cytosol due to transport form
the extracellular space, or due to degradation of certain P-containing biomolecules. Nucleoside
triphosphates (e.g. ATP), nucleic acid (DNA or RNA), and membrane lipids (e.g. phospholipids)
are all biomolecular sources of P. The mentioned sources of P are of course vital for all life on Earth,
meaning that there must be tight regulations for P uptake and management within living cells.
Some of the proteins that regulate this are Pi transporters, phosphatases, and phospholipases.
The former are proteins that transport Pi across membranes. Phosphatases and phospholipases
are both groups of enzymes that make P more accessible for the synthesis of new compounds, by
cleaving covalent bonds in P-containing biomolecules

# Experiment description
P. tricornutum cells are grown in two different growth media, one with surplus of phosphate,
(phosphate replete) and one with reduced phosphate levels (phosphate deplete).
• RNA is isolated from three biological replicas from each, replete or deplete, after 48 hours growth.
• As the cells divide the concentration of phosphate in the growth medium is reduced, and the Pdepleted cells starts to have problems.
• Too little phosphate will reduce the production of nucleotides which in turn affects DNA synthesis
and reduce the rate of cell division.
• The cells start processes that increase uptake of phosphate (phosphate transporters) and
increase scavenging of phosphate from external and internal sources (various phosphatases).
The lipid profile is also affected, and levels of phospholipids are reduced.

# Program description
The R-script will use the functions included in edgeR package.
The R-script will perform several operations and analyses:
• The files containing the count table (read counts) and sample information are read into R.
• The count table will be filtered so that genes with very few reads are excluded from the
analysis.
• The data will be analyzed using statistical methods based on generalized linear models
(GLM).
• It will perform likelihood ratio tests to find differentially expressed genes
• It will produce a plot which gives an overview of how many genes are differentially expressed
• It will produce a plot which shows how much biological variation you have in the biological
replicas.

