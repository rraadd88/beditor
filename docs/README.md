beditor
=======

A computational workflow for designing libraries of guide RNAs for
CRISPR base editing

Installation
------------

1.  Create a virtual environment.

``` {.sourceCode .text}
cd beditor
conda create -f environment.yml
```

2.  Activate the virtual environment.

``` {.sourceCode .text}
source activate beditor
```

3.  Run the analysis.

``` {.sourceCode .text}
beditor --cfg configuration.yml
```

Configuration file (configuration.yml)
--------------------------------------

It contains the all the analysis specific parameters.

Template:
<https://github.com/rraadd88/test_beditor/blob/master/common/configuration.yml>

Input format
------------

``` {.sourceCode .text}
`mutation_format` opted in configuration.yml file and corresponding columns needed in input: 
`nucleotide` : `['genome coordinate','nucleotide mutation']`.
`aminoacid`  : `['transcript: id','aminoacid: position','amino acid mutation']`.
```

Output format
-------------

``` {.sourceCode .text}
`mutation_format` opted in configuration.yml file and corresponding columns needed in input: 

`nucleotide` :  ['genome coordinate','nucleotide wild-type','nucleotide mutation',]
`nucleotide` : ['transcript: id','aminoacid: wild-type','aminoacid: position','amino acid mutation','codon: wild-type','guide: id','guide+PAM sequence','beditor score','alternate alignments count','CFD score']
```

Format of \`guide: id\`:

``` {.sourceCode .text}
{genomic locus}|{position}|({strategy})
where,
strategy= {base editor};{strand};@{distance of mutation from PAM};{PAM};{codon wild-type}:{codon mutation};{amino acid wild-type}:{amino acid mutation};
```

A directory by the basename of configuration file (eg. directory called
'human' if configuration file is 'human.yml') would be created in the
same folder where configuration file is located. It is referred to as
'project directory'.

Inside a project directory there would be following folders named by
corresponding steps of analysis.

``` {.sourceCode .text}
1. `01_sequences/`
Stores the output of step #1. Extracting sequences flanking mutation site.
2. `02_mutagenesis/`
Stores the output of step #2. Estimating the editable mutations based on base editors chosen.
3. `03_guides/`
Stores the output of step #3. Designed guides.
4. `04_offtargets/`
Stores the output of step #4. Offtarget effects.
5. `05_outputs/`
Stores combined output and visualizations.

Also,
- `00_input/`
Store the input files.
- `chunks/`
If parallel processing is used, this folder would store individual parts (chunks) of the analysis. 
```

How to install new base editor or PAM
-------------------------------------

To list supported pams and editors

``` {.sourceCode .text}
beditor --list pams
beditor --list editors
```

These lists are located in beditor/data folder and can be modified.

Working with non-ensembl genome
-------------------------------

<https://github.com/openvax/pyensembl#non-ensembl-data>

How to analyze test dtasets
---------------------------

``` {.sourceCode .text}
# make the input files with mock data
git clone https://github.com/rraadd88/test_beditor.git
source activate beditor;cd test_beditor;python make_datasets.py
# command to analyze mock input data of a species
source activate beditor;cd dataset_{species name};beditor --cfg mutation_format_nucleotide_mutation_mutations_for.yml
```

API
---
