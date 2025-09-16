Datasets Included
=================

This page provides an overview of the datasets included in CoderData version 2.1.0.

Figshare record: https://api.figshare.com/v2/articles/28823159
Version: 2.1.0

---------------------------
Dataset Overview
---------------------------
.. csv-table:: Datasets and Modalities
   :header: "Dataset", "References", "Sample", "Transcriptomics", "Proteomics", "Mutations", "Copy Number", "Drug", "Drug Descriptor", "Experiments"
   :widths: 12, 10, 6, 12, 12, 12, 12, 8, 15, 12

   "BeatAML", "[1]_, [2]_", "X", "X", "X", "X", "", "X", "X", "X"
   "BladderPDO", "[3]_", "X", "X", "", "X", "X", "X", "X", "X"
   "CCLE", "[4]_", "X", "X", "X", "X", "X", "X", "X", "X"
   "CPTAC", "[5]_", "X", "X", "X", "X", "X", "", "", ""
   "CTRPv2", "[6]_, [7]_, [8]_", "X", "X", "", "X", "X", "X", "X", "X"
   "FIMM", "[9]_, [10]_", "X", "X", "", "", "", "X", "X", "X"
   "HCMI", "[11]_", "X", "X", "", "X", "X", "", "", ""
   "MPNST", "[12]_", "X", "X", "X", "X", "X", "X", "X", "X"
   "NCI60", "[13]_", "X", "X", "X", "X", "", "X", "X", "X"
   "Pancreatic PDO", "[14]_", "X", "X", "", "X", "X", "X", "X", "X"
   "PRISM", "[15]_, [16]_", "X", "X", "", "", "", "X", "X", "X"
   "Sarcoma PDO", "[17]_", "X", "X", "", "X", "", "X", "X", "X"
   "CRC PDO", "[18]_", "X", "X", "", "X", "X", "X", "X", ""
   "Liver PDO", "[19]_", "X", "X", "", "X", "X", "X", "X", ""
   "Novartis PDX", "[20]_", "X", "X", "", "X", "X", "X", "X", ""
   "gCSI", "[21]_, [22]_", "X", "X", "X", "X", "X", "X", "X", ""
   "GDSC v1", "[23]_, [24]_, [25]_", "X", "X", "X", "X", "X", "X", "X", ""
   "GDSC v2", "[23]_, [24]_, [25]_", "X", "X", "X", "X", "X", "X", "X", ""

The table above lists the datasets included in CoderData version 2.1.0, along with references to their original publications and the types of data available for each dataset. An "X" indicates the presence of a particular data type for the corresponding dataset.


---------------------------
Dataset Summary Statistics
---------------------------
The following table summarizes key statistics for each dataset, including the number of samples, drugs, and various combinations of sample-drug pairs with different molecular data types.

    .. csv-table:: Dataset Summary Statistics
       :file: _static/dataset_summary_statistics.csv
       :header-rows: 0


---------------------------------
Drug Curve Metrics Collected
---------------------------------
The following table summarizes the number of drugs associated with each dose-response metric across the datasets.

    .. csv-table:: Drug Curve Metrics Summary
       :file: _static/dataset_curve_metric_summary.csv
       :header-rows: 0



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
