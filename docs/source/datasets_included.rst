Datasets Included
=================

This page provides an overview of the datasets included in CoderData version 2.2.0. This package collects 18 diverse sets of paired molecular datasets with corresponding drug sensitivity data. All data here is reprocessed and standardized so it can be easily used as a benchmark dataset for machine learning models.

Figshare record: https://api.figshare.com/v2/articles/28823159
Version: 2.2.0

---------------------------
Dataset Overview
---------------------------
.. csv-table:: Datasets and Modalities
   :header: "Dataset", "References", "Sample", "Drug", "Drug Descriptor", "Experiments", "Transcriptomics", "Proteomics", "Mutations", "Copy Number"
   :widths: 14, 12, 6, 8, 15, 12, 12, 12, 12, 12

   "BeatAML", "[1]_, [2]_", "1022", "164", "X", "X", "X", "X", "X", ""
   "Bladder", "[3]_", "134", "50", "X", "X", "X", "", "X", "X"
   "CCLE", "[4]_", "502", "24", "X", "X", "X", "X", "X", "X"
   "Colorectal ", "[18]_", "61", "10", "X", "", "X", "", "X", "X"
   "CPTAC", "[5]_", "1139", "", "", "", "X", "X", "X", "X"
   "CTRPv2", "[6]_, [7]_, [8]_", "846", "459", "X", "X", "X", "", "X", "X"
   "FIMM", "[9]_, [10]_", "52", "52", "X", "X", "X", "", "", ""
   "GDSC v1", "[23]_, [24]_, [25]_", "984", "294", "X", "", "X", "X", "X", "X"
   "GDSC v2", "[23]_, [24]_, [25]_", "806", "171", "X", "", "X", "X", "X", "X"
   "gCSI", "[21]_, [22]_", "569", "X", "X", "", "X", "X", "X", "X"
   "HCMI", "[11]_", "886", "", "", "", "X", "", "X", "X"
   "Liver", "[19]_", "62", "76", "X", "", "X", "", "X", "X"
   "MPNST", "[12]_", "50", "30", "X", "X", "X", "X", "X", "X"
   "NCI60", "[13]_", "83", "55157", "X", "X", "X", "X", "X", ""
   "Novartis", "[20]_", "386", "25", "X", "", "X", "", "X", "X"
   "Pancreatic", "[14]_", "70", "25", "X", "X", "X", "", "X", "X"
   "PRISM", "[15]_, [16]_", "478", "1419", "X", "X", "X", "", "", ""
   "Sarcoma", "[17]_", "36", "34", "X", "X", "X", "", "X", ""


The table above lists the datasets included in CoderData version 2.2.0, along with references to their original publications, counts of samples and drugs, and the types of data available for each dataset.

CoderData includes the following data:

- Sample - cell lines, patient-derived samples, or patient-derived organoids
- Drug - compounds tested for sensitivity
- Drug Descriptor - molecular descriptors for each drug (computed using RDKit)
- Experiments - dose-response experiments (various metrics such as AUC, IC50, etc.)
- Transcriptomics - gene expression (in transcripts per million, TPM)
- Proteomics - protein expression (in log2 ratio to reference)
- Mutations - gene mutations (variant calls)
- Copy Number - gene copy number variations (number of copies of each gene, 2 being diploid)

An "X" indicates the presence of a particular data type for the corresponding dataset.


---------------------------
Dataset Summary Statistics
---------------------------
The following table summarizes combination counts for each dataset. This includes the number of experimental sample-drug pairs, with different molecular data types. Each column represents the number of unique combinations of samples and drugs with the specified molecular data types available. For example, the "Sample-Drug-Transcriptomics-Mutations" column indicates the number of unique sample-drug pairs that have both transcriptomics and mutation data available.

    .. csv-table:: Dataset Summary Statistics
       :file: _static/dataset_summary_statistics.csv
       :header-rows: 0


---------------------------------
Drug Curve Metrics Collected
---------------------------------
The following table summarizes the number of drugs associated with each dose-response metric across the datasets.

    .. csv-table:: Drug Curve Metrics Summary
       :file: _static/dataset_curve_metrics_wide.csv
       :header-rows: 0

Types of dose-response metrics collected include:

- AAC - Area above the response curve; the complement value of AUC.
- ABC - Area between curves, the difference between the AUC of the control and the treated cells.
- AUC - Area under the fitted hill slope curve across all doses present. Lower AUC signifies lower levels of growth.
- DSS - A multiparametric dose response value that takes into account control and treated cells.
- fit_auc - Area under the fitted hill slope curve across the common interval of −log10[M], where the molar concentration ranges from 10⁻⁴ to 10⁻¹⁰.
- fit_ec50 - The fitted curve prediction of the −log10M concentration at which 50% of the maximal effect is observed.
- fit_ec50se - Standard error of the Fit_EC50 estimate.
- fit_einf - The fraction of cells that are unaffected even at an infinite dose concentration. Calculated as the lower asymptote of the hill slope function.
- fit_hs - The estimated hill slope binding cooperativity, calculated as the slope of the sigmoidal hill curve.
- fit_ic50 - The fitted curve prediction of the −log10M concentration required to reduce tumor growth by 50%.
- fit_r2 - Coefficient of determination between observed growth and the fitted hill slope curve, indicating goodness of fit.
- lmm - The resulting “time and treatment interaction” in a linear mixed model with fixed effects as time and treatment and patient as a random effect. Indicates how much the treatment changes the slope of log(volume) over time compared to the control.
- mRESCIST - Disease status classified into PD (progressive disease), SD (stable disease), PR (partial response), and CR (complete response), based on percent volume change and cumulative average response.
- published_auc - Published Area Under the Curve
- TG - Tumor growth inhibition between the control and treatment time-volume curves.



---------------------------
References
---------------------------

.. [1] Bottomly D, Long N, Schultz AR, et al. *Integrative analysis of drug response and clinical outcome in acute myeloid leukemia.* Cancer Cell. 2022;40(8):850-864.e9. doi:`10.1016/j.ccell.2022.07.002 <https://doi.org/10.1016/j.ccell.2022.07.002>`_
.. [2] Pino JC, Posso C, Joshi SK, et al. *Mapping the proteogenomic landscape enables prediction of drug response in acute myeloid leukemia.* Cell Rep Med. 2024;5(1):101359. doi:`10.1016/j.xcrm.2023.101359 <https://doi.org/10.1016/j.xcrm.2023.101359>`_
.. [3] Lee SH, Hu W, Matulay JT, et al. *Tumor Evolution and Drug Response in Patient-Derived Organoid Models of Bladder Cancer.* Cell. 2018;173(2):515-528.e17. doi:`10.1016/j.cell.2018.03.017 <https://doi.org/10.1016/j.cell.2018.03.017>`_
.. [4] Barretina J, Caponigro G, Stransky N, et al. *The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity.* Nature. 2012;483(7391):603-607. doi:`10.1038/nature11003 <https://doi.org/10.1038/nature11003>`_
.. [5] Lindgren CM, Adams DW, Kimball B, et al. *Simplified and Unified Access to Cancer Proteogenomic Data.* J Proteome Res. 2021;20(4):1902-1910. doi:`10.1021/acs.jproteome.0c00919 <https://doi.org/10.1021/acs.jproteome.0c00919>`_
.. [6] Rees MG, Seashore-Ludlow B, Cheah JH, et al. *Correlating chemical sensitivity and basal gene expression reveals mechanism of action.* Nat Chem Biol. 2016;12(2):109-116. doi:`10.1038/nchembio.1986 <https://doi.org/10.1038/nchembio.1986>`_
.. [7] Seashore-Ludlow B, Rees MG, Cheah JH, et al. *Harnessing Connectivity in a Large-Scale Small-Molecule Sensitivity Dataset.* Cancer Discov. 2015;5(11):1210-1223. doi:`10.1158/2159-8290.CD-15-0235 <https://doi.org/10.1158/2159-8290.CD-15-0235>`_
.. [8] Basu A, Bodycombe NE, Cheah JH, et al. *An interactive resource to identify cancer genetic and lineage dependencies targeted by small molecules.* Cell. 2013;154(5):1151-1161. doi:`10.1016/j.cell.2013.08.003 <https://doi.org/10.1016/j.cell.2013.08.003>`_
.. [9] Mpindi JP, Yadav B, Östling P, et al. *Consistency in drug response profiling.* Nature. 2016;540(7631):E5-E6. doi:`10.1038/nature20171 <https://doi.org/10.1038/nature20171>`_
.. [10] Pemovska T, Kontro M, Yadav B, et al. *Individualized systems medicine strategy to tailor treatments for patients with chemorefractory acute myeloid leukemia.* Cancer Discov. 2013;3(12):1416-1429. doi:`10.1159/2159-8290.CD-13-0350 <https://doi.org/10.1158/2159-8290.CD-13-0350>`_
.. [11] Human Cancer Models Initiative (HCMI). dbGaP accession phs001486. `https://cancer.gov/ccg/research/functional-genomics/hcmi <https://cancer.gov/ccg/research/functional-genomics/hcmi>`_
.. [12] Dehner C, Moon CI, Zhang X, et al. *Chromosome 8 gain is associated with high-grade transformation in MPNST.* JCI Insight. 2021;6(6):e146351. doi:`10.1172/jci.insight.146351 <https://doi.org/10.1172/jci.insight.146351>`_
.. [13] Shoemaker RH. *The NCI60 human tumour cell line anticancer drug screen.* Nat Rev Cancer. 2006;6(10):813-823. doi:`10.1038/nrc1951 <https://doi.org/10.1038/nrc1951>`_
.. [14] Tiriac H, Belleau P, Engle DD, et al. *Organoid Profiling Identifies Common Responders to Chemotherapy in Pancreatic Cancer.* Cancer Discov. 2018;8(9):1112-1129. doi:`10.1158/2159-8290.CD-18-0349 <https://doi.org/10.1158/2159-8290.CD-18-0349>`_
.. [15] Corsello SM, Nagari RT, Spangler RD, et al. *Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling.* Nat Cancer. 2020;1(2):235-248. doi:`10.1038/s43018-019-0018-6 <https://doi.org/10.1038/s43018-019-0018-6>`_
.. [16] Yu C, Mannan AM, Yvone GM, et al. *High-throughput identification of genotype-specific cancer vulnerabilities in mixtures of barcoded tumor cell lines.* Nat Biotechnol. 2016;34(4):419-423. doi:`10.1038/nbt.3460 <https://doi.org/10.1038/nbt.3460>`_
.. [17] Al Shihabi A, Tebon PJ, Nguyen HTL, et al. *The landscape of drug sensitivity and resistance in sarcoma.* Cell Stem Cell. 2024;31(10):1524-1542.e4. doi:`10.1016/j.stem.2024.08.010 <https://doi.org/10.1016/j.stem.2024.08.010>`_
.. [18] van de Wetering M, Francies HE, Francis JM, et al. *Prospective derivation of a living organoid biobank of colorectal cancer patients.* Cell. 2015;161(4):933-945. doi:`10.1016/j.cell.2015.03.053 <https://doi.org/10.1016/j.cell.2015.03.053>`_
.. [19] Ji S, Feng L, Fu Z, et al. *Pharmaco-proteogenomic characterization of liver cancer organoids for precision oncology.* Sci Transl Med. 2023;15(706):eadg3358. doi:`10.1126/scitranslmed.adg3358 <https://doi.org/10.1126/scitranslmed.adg3358>`_
.. [20] Gao H, Korn JM, Ferretti S, et al. *High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response.* Nat Med. 2015;21(11):1318–1325. doi:`10.1038/nm.3954 <https://doi.org/10.1038/nm.3954>`_
.. [21] Haverty PM, Lin E, Tan J, et al. *Reproducible pharmacogenomic profiling of cancer cell line panels.* Nature. 2016;533(7603):333–337. doi:`10.1038/nature17987 <https://doi.org/10.1038/nature17987>`_
.. [22] Klijn C, Durinck S, Stawiski EW, et al. *A comprehensive transcriptional portrait of human cancer cell lines.* Nat Biotechnol. 2015;33(3):306–312. doi:`10.1038/nbt.3080 <https://doi.org/10.1038/nbt.3080>`_
.. [23] Garnett MJ, Edelman EJ, Heidorn SJ, et al. *Systematic identification of genomic markers of drug sensitivity in cancer cells.* Nature. 2012;483(7391):570–575. doi:`10.1038/nature11005 <https://doi.org/10.1038/nature11005>`_
.. [24] Iorio F, Knijnenburg TA, Vis DJ, et al. *A Landscape of Pharmacogenomic Interactions in Cancer.* Cell. 2016;166(3):740–754. doi:`10.1016/j.cell.2016.06.017 <https://doi.org/10.1016/j.cell.2016.06.017>`_
.. [25] Yang W, Soares J, Greninger P, et al. *Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells.* Nucleic Acids Res. 2013;41(Database issue):D955–D961. doi:`10.1093/nar/gks1111 <https://doi.org/10.1093/nar/gks1111>`_
