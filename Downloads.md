**NOTE** - parts of code comments need cleaning. Please feel free to use the code and make amends as currently i do not have Matlab license.

# Files #

---


  * twoHoldOutExp (m file)
  * readCustomFile (m file)
  * dataStorage (m file)
  * generateInteraction (m file)
  * generateGeneCpd (m file)
  * geneTRCMPLXstats (m file)
  * plotAUC (m file)
  * geneExpression (mat file)
  * probabilities (mat file)
  * aucANDpredictions sample TRCMPLX.mat (mat file)

# Description of above files #

---


  * _twoHoldOutExp_

> conducts the two hold out experiment and is the main file/function to start with. Takes _eviDence_ ('ge' for gene expression, 'me' for methylation) and _model_ ('t1' for PBK+EI, 't2' for PBK and 'p1' for NB+MPBK). Here, PBK+EI - Prior Biological Knowledge and Epigenetic Information, PBK - Prior Biological Knowledge, NB+MPBK - Naive Bayes with Minimanl Prior Biological Knowledge

  * _readCustomFile_

> reads the gene expression from the _geneExpression_ mat file and returns a set of unique genes in _uniqueGenes_, the gene expression matrix in _expressionMatrix_, the total number of genes involved _noGenes_, the total number of samples involved _noSamples_, the ground truth labels of the samples in _groundTruthLabels_ and _transGroundTruthLabels_

  * _dataStorage_

> take the type of model involved stored in _model_ ('t1'/'t2'/'p1') and stores the required probability values in workspace as well as _probabilities_ mat file.

  * _generateInteraction_

> take the unique genes stored in _uniqueGenes_ as well as the _model_ ('t1'/'t2'/'p1') and generates the interaction (i.e causal arcs) among the nodes in _interaction_ and the list of names of the nodes in _nodeNames_. These will be used in the generation of the directed acyclic graph of the Bayesian Network.

  * _generateGenecpd_

> takes in the expression data of a particular gene for training in _vecTraining_, the labels pertaining to the particular gene expression data in _labelTraining_, the name of the gene node in _nodeName_, the parents of the gene node in _parent_, the probability tables of the parents of the gene node in _parent\_cpd_ and finally the model type in _model_ ('t1'/'t2'/'p1') to return the conditional probability table for the gene node under consideration.

  * _geneTRCMPLXstats_

> This file generates the tables which shows how TRCMPLX behaves as the evidences of genes vary in both normal and tumorous cases and are depicted in the results sections of Sinha:2014. Interpretations of the results can be studied in more depth from Sinha:2014. Once the results have been saved in Results.mat file, one can rename the file based on the model and the evidence arguments used in function twoHoldOutExp. An example of one such file is Results-T1-GE-pforTRCMPLX-90per.mat which is loaded using the load command and further processed using the script in the m file titled geneTRCMPLXstats.

  * _plotAUC_

> Note that to generate the ROC graphs and their respective AUC values for different models with varying effect of TRCMPLX on different genes (ETGN in Sinha:2014), the results in variable s X and Y (of twoHoldOutExp) are stored in different variables and clumped together in a mat file titled aucANDpredictions\_sample\_TRCMPLX.mat. This has to be done manually for each model and every setting of ETGN. For example, using model t1 and ETGN of 60%, the false positive rate in X is stored as xT1\_50 and the true positive rate in Y is stored as yT1\_50, in the above mentioned .mat file. Finally, the script in the m file titled plotAUC is used to manipulate the aforementioned transformed variables and generate the ROC curves in the results section of Sinha:2014.

# Files needed to run #

---


> Before running the scripts
  1. one must have Matlab runing on the machine
  1. download and install the Bayesian Network Toolboox by Murphy on [[BNT-link](https://code.google.com/p/bnt/)]

# Bibliography #

---


In case of use of the code please cite the articles in (1) and (2).

  1. _Integration of prior biological knowledge and epigenetic information enhances the prediction accuracy of the Bayesian Wnt pathway_; Shriprakash Sinha; Integr. Biol., 2014, 6(11), 1034-1048, The Royal Society of Chemistry [[RSC-link](http://dx.doi.org/10.1039/C4IB00124A)]

```
@Article{C4IB00124A,
author ="Sinha, Shriprakash",
title  ="Integration of prior biological knowledge and epigenetic information enhances the prediction accuracy of the Bayesian Wnt pathway",
journal  ="Integr. Biol.",
year  ="2014",
volume  ="6",
issue  ="11",
pages  ="1034-1048",
publisher  ="The Royal Society of Chemistry",
doi  ="10.1039/C4IB00124A",
url  ="http://dx.doi.org/10.1039/C4IB00124A",
}
```

  1. _A pedagogical walkthrough of computational modeling and simulation of Wnt signaling pathway using static causal models in Matlab_; Shriprakash Sinha; Unpublished Manuscript, 2014, 1-25, doi : 10.1101/011064 [[Biorxiv-link](http://biorxiv.org/content/early/2015/01/29/011064)]

  1. Please cite paper by Jiang et.al. 2008 on _DACT3 is an epigenetic regulator of Wnt/beta-catenin signaling in colorectal cancer and is a therapeutic target of histone modifications_ for gene expression dataset.