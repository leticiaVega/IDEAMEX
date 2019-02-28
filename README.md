# IDEAMEX

IDEAMEX (Integrative Differential Expression Analysis for Multiple EXperiments) is a web server where
users can run simultaneously, the best Bioconductor (Huber W 2015) packages for RNA-seq differential
expression analysis. The web server also integrates the results from each package, so users can select between
the intersection or union of all results. The sole input for the IDEAMEX pipeline, is a raw count table in
simple text format, for as many desired replicates and conditions, allowing the user to select which conditions
will be compared, according to the biological design behind the experiment. The process consists of three
main steps 1) Data Analysis: that allows a preliminary analysis for data quality control based on the data
distribution per sample, using different types of graphs; 2) Differential expression: performs the differential
expression analysis using the bioconductor packages: edgeR (Robinson MD 2010), limma - voom (Ritchie
ME 2015), DESeq2 (Love MI 2014) and NOISeq (D. J. Tarazona S Garcia-Alcalde F 2011), and generates
reports for each method; 3) Result integration: the integrated results are reported using Venn diagrams,
heatmaps, correlograms and text lists where differentially expressed genes are reported, according to the
cutoff lines defined by the user. Our server allows an easy and friendly visualization of results, providing an
easy interaction during the analysis process, as well as error tracking and debugging by providing output log
files.
