# SAGAconf User Manual

## Introduction

SAGAconf is a software that assigns calibrated confidence scores to chromatin state annotations, addressing the problem of reproducibility of SAGA annotations. It uses an integrative approach to derive a reproducibility score (r-value) at each genomic position, allowing for the identification of a confident and reliable subset from genome annotation while removing irreproducible predictions. Thus, SAGAconf allows a researcher to select only the reliable predictions from a chromatin annotation for use in downstream analyses. SAGAconf is independent of SAGA methodology, taking as input only posterior probability matrix.

## Input Data

SAGAconf reproducibility analysis pipeline gets as input two sets of posterior probability matrices (base and verification) each with **`G`** rows and **`K+3`** columns. Here, **`G`** denotes the genome size that is the number of rows and **`K`** corresponds to the number of chromatin states identified by SAGA method with 3 additional columns representing genomic coordinates of each bin (**`chr`**, **`start`**, **`end`**). These posterior files can be in BED or CSV format.

---

---

The **`SAGAconf_parser.py`** script is a command-line tool that allows you to parse posterior files generated by the SAGA model (ChromHMM and Segway). The script takes several arguments, some of which are required and others are optional.

## Required Arguments

- **`posteriordir`**: This argument specifies the directory containing all posterior files. You must provide the path to this directory when running the script.
- **`resolution`**: This argument specifies the resolution of the SAGA model in base pairs (bp). You must provide this value when running the script.
- **`savedir`**: This argument specifies the directory where the parsed posterior file will be saved. You must provide the path to this directory when running the script.
- **`--saga`**: This argument specifies which SAGA model to use. You can choose between **`segway`** and **`chmm`**.

## Optional Arguments

- **`-h`** or **`--help`**: This argument displays a help message and exits the script. You can use this argument if you need more information about how to use the script.
- **`--out-format`**: This argument specifies the format for saving parsed posteriors. You can choose between **`bed`** and **`csv`** formats. If you do not provide this argument, the default format **`bed`**  will be used.

## Example Usage

Here is an example of how you might run the **`SAGAconf_parser.py`** script:

`python SAGAconf_parser.py --out_format bed --saga chmm /path/to/posteriordir 200 /path/to/savedir`

In this example, we are specifying that we want to use the **`bed`** format for saving parsed posteriors, and we want to use the **`chmm`** SAGA model. We are also providing the required arguments for the **`posteriordir`**, **`resolution`**, and **`savedir`**.

---

---

Given the parserd_posterior files, SAGAconf.py can provide following analysis:

1. ratio of naive overlap between chromatin states of base and verification annotations
2. Granularity of chromatin states vs. overlap
3. Spatial misalignment of chromatin states vs. overlap
4. Posterior probability vs. overlap
5. Post-clustering of chromatin states by merging into <K chromatin states
6. r-Values per genomic position
7. Robust chromatin states (confident segments) 

---

---

The **`SAGAconf.py`** script is a command-line tool that allows you to compare base and verification annotations generated by the SAGA model. The script takes several arguments, some of which are required and others are optional.

## Required Arguments

- **`base`**: This argument specifies the path of the parsed posterior CSV or BED file of the base replicate annotation. You must provide the path to this file when running the script.
- **`verif`**: This argument specifies the path of the parsed posterior CSV or BED file of the verification replicate annotation. You must provide the path to this file when running the script.
- **`savedir`**: This argument specifies the directory where the SAGAconf results will be saved. You must provide the path to this directory when running the script.

## Optional Arguments

- **`-h`** or **`--help`**: This argument displays a help message and exits the script. You can use this argument if you need more information about how to use the script.
- **`-v`** or **`--verbosity`**: This argument increases output verbosity. You can use this argument if you want more detailed output from the script.
- **`-bm`** or **`--base_mnemonics`**: This argument specifies a file to use as mnemonics (biological label interpretations) for the base replicate. If you provide this argument, you must also provide the path to the text file.
- **`-vm`** or **`--verif_mnemonics`**: This argument specifies a file to use as mnemonics (biological label interpretations) for the verification replicate. If you provide this argument, you must also provide the path to the text file.
- **`-s`** or **`--subset`**: This argument specifies that SAGA should be run on just one chromosome (chr21). If you provide this argument, only one chromosome (chr21) will be analyzed. Note that, if your annotations do not include chr21, you should not use this argument.
- **`-w`** or **`--windowsize`**: This argument specifies the window size in base pairs (bp) to account for around each genomic bin. The default value is 1000 bp.
- **`-to`** or **`--iou_threshold`**: This argument specifies the threshold on the Intersection over Union (IoU) of overlap for considering a pair of labels as corresponding. You can provide a value for this threshold when running the script. The default value is 0.75.
- **`-tr`** or **`--repr_threshold`**: This argument specifies the threshold on the reproducibility score for considering a segment as reproduced. You can provide a value for this threshold when running the script. In the original paper, this argument is referred to as **`alpha`.** The default value is 0.8.
- **`-k`** or **`--merge_clusters`**: This argument specifies a k value to merge base annotation chromatin states until k states. You can provide a value for k when running the script.
- **`-q`** or **`--quick`**: This argument specifies that only a subset of essential analysis should be performed for a quick report. If you provide this argument, only essential analysis will be performed.

## Example Usage

Here is an example of how you might run the **`SAGAconf.py`** script:

`python SAGAconf/SAGAconf.py -v -bm /path/to/base_mnemonics.txt -vm /path/to/verif_mnemonics.txt -s -w 2000 -to 0.5 -tr 0.8 -k 5 -q /path/to/base/parsed_posterior.bed /path/to/verif/parsed_posterior.bed /path/to/savedir`

In this example, we are specifying that we want increased output verbosity (**`-v`**), and we are providing paths to text files containing mnemonics for both base and verification replicates (**`-bm`** and **`-vm`**). We are also specifying that we want to run SAGA on just chr21 chromosome (**`-s`**), and we are providing values for window size (**`-w`**), IoU threshold (**`-to`**), reproducibility threshold (**`-tr`**), and k value (**`-k`**). Finally, we are specifying that only essential analysis should be performed (**`-q`**) and providing paths to base and verification replicate files and save directory.

## Mnemonics file

A mnemonics file is a tab-separated text file with two columns. The first column contains the integer label of the original chromatin states, and the second column contains the biological label assignments. Here is an example of a mnemonics file:


```
old	new
1	Enhancer_low
2	Enhancer
3	Promoter_flanking
4	Promoter
5	Promoter
6	Promoter
7	Enhancer
8	Enhancer
9	Enhancer
10	Enhancer
```


To generate your own mnemonics file, you can create a new text file and enter the integer labels and biological label assignments in the appropriate columns. Make sure to separate the columns with a tab character and to include a header row with the column names **`old`** and **`new`**.

---

---

## Real Use Case Example

In this section, we will go through a real use case of SAGAconf starting from obtaining chromatin state annotations from ChromHMM, followed by parsing the posteriors and running SAGAconf. Note that, in the following example, everything is happening within the main `SAGAconf/` directory.

1. First, download ChromHMM and unzip it by running the following commands:
    
    ```
    curl -o ChromHMM.zip http://compbio.mit.edu/ChromHMM/ChromHMM.zip
    unzip ChromHMM.zip
    ```
    
2. Next, generate test annotations from ChromHMM by running the following command:
    
    `java -mx1600M -jar ChromHMM/ChromHMM.jar LearnModel -printposterior ChromHMM/SAMPLEDATA_HG18 ChromHMM/OUTPUTSAMPLE 10 hg18`
    
    This command runs the `LearnModel` command of the `ChromHMM` program with the `-printposterior` option. This command takes a set of binarized data files, learns chromatin state models, and produces chromatin state annotations with posterior files where `k=10`. In the ChromHMM sample run, it concatenates the data from `K562_chr11` and `GM12878_chr11`. For demonstration purposes, we will treat these two as base and verification replicates.
    
3. Next, we manually move the posterior files for `GM12878_chr11` and `K562_chr11` into separate directories for base and verification replicates using the following command:

    ```
    mkdir -p ChromHMM/OUTPUTSAMPLE/base ChromHMM/OUTPUTSAMPLE/verif && \
    mv ChromHMM/OUTPUTSAMPLE/POSTERIOR/GM12878_10_chr11_posterior.txt ChromHMM/OUTPUTSAMPLE/base && \
    mv ChromHMM/OUTPUTSAMPLE/POSTERIOR/K562_10_chr11_posterior.txt ChromHMM/OUTPUTSAMPLE/verif && \
    rm -r ChromHMM/OUTPUTSAMPLE/POSTERIOR
    ```
    
    This command creates two directories, `ChromHMM/OUTPUTSAMPLE/base` and `ChromHMM/OUTPUTSAMPLE/verif`, moves the posterior files for `GM12878_chr11` and `K562_chr11` into these directories, respectively, and removes the original `POSTERIOR` directory.
    
4. Then, we use the following commands to parse posteriors from `ChromHMM/OUTPUTSAMPLE/base` and `ChromHMM/OUTPUTSAMPLE/verif` into the standard format required by SAGAconf:
    
    ```
    python SAGAconf_parser.py --saga chmm ChromHMM/OUTPUTSAMPLE/base 200 ChromHMM/OUTPUTSAMPLE/base

    python SAGAconf_parser.py --saga chmm ChromHMM/OUTPUTSAMPLE/verif 200 ChromHMM/OUTPUTSAMPLE/verif
    ```
    
    Here, using `--saga chmm`, we are specifying that the annotations are generated by ChromHMM and we are specifying resolution to be 200bp. Ultimately, we will have `parsed_posterior.bed` files in both `ChromHMM/OUTPUTSAMPLE/base` and `ChromHMM/OUTPUTSAMPLE/verif`.
    
5. Lastly, we run SAGAconf to obtain a full reproducibility report with default parameters by running the following command:
    
    `python SAGAconf/SAGAconf.py ChromHMM/OUTPUTSAMPLE/base/parsed_posterior.bed ChromHMM/OUTPUTSAMPLE/verif/parsed_posterior.bed ChromHMM/OUTPUTSAMPLE/sagaconf_base`