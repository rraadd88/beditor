# `beditor`

Python package to streamline designing of guides. 

The package (called `beditor`) can be installed in python 3.6 environment by following commands.

```
cd beditor
pip install -e .
```

This command will run an analysis:

```
beditor <path to a json configuration file>
```

Here, json configuration file contains path to a csv file with transcript ids and mutation position (residue number) and other parameters. (See `beditor --help` for more options. )

beditor/test folder contains a demo configuration file and mock data that can be run by following commands:

```
cd beditor/test
beditor human_test.json
```

Note that when running for the first time, pyensembl would download genome data (91th release,~400Mb). It would take some time.
Also, I have not connected the modules for yeast yet, so currently human data can analysed using this package.


## Tests

### For Human sequences

	(01pppi)$ beditor human_test.json
	[2018-04-23 08:57:34,649] INFO	from pipeline.py in main(..):42: start
	[2018-04-23 08:57:45,674] INFO	from sequence_data.py in _load_or_create_fasta_dictionary_pickle(..):118: Loaded sequence dictionary from /home/anonymous/.cache/pyensembl/GRCh38/ensembl91/Homo_sapiens.GRCh38.cdna.all.fa.gz.pickle
	[2018-04-23 08:57:46,039] INFO	from sequence_data.py in _load_or_create_fasta_dictionary_pickle(..):118: Loaded sequence dictionary from /home/anonymous/.cache/pyensembl/GRCh38/ensembl91/Homo_sapiens.GRCh38.ncrna.fa.gz.pickle
	[2018-04-23 08:57:46,779] INFO	from sequence_data.py in _load_or_create_fasta_dictionary_pickle(..):118: Loaded sequence dictionary from /home/anonymous/.cache/pyensembl/GRCh38/ensembl91/Homo_sapiens.GRCh38.pep.all.fa.gz.pickle
	[2018-04-23 08:57:50,545] INFO	from get_seq.py in din2dseq(..):124: Counts of amino acids to mutate:
	[2018-04-23 08:57:50,547] INFO	from get_seq.py in din2dseq(..):125: S    34
	T    33
	Y    32
	Name: aminoacid: wild-type, dtype: int64
	Possible 1 nucleotide mutations of the phospho-sites:
	           amino acid mutation      method codon codon mutation
	amino acid                                                     
	T                            A         ABE   ACA            GCA
	T                            A         ABE   ACC            GCC
	T                            A         ABE   ACG            GCG
	T                            A         ABE   ACT            GCT
	S                            G         ABE   AGC            GGC
	S                            G         ABE   AGT            GGT
	Y                            H         ABE   TAC            CAC
	Y                            C         ABE   TAC            TGC
	Y                            H         ABE   TAT            CAT
	Y                            C         ABE   TAT            TGT
	S                            P         ABE   TCA            CCA
	S                            P         ABE   TCC            CCC
	S                            P         ABE   TCG            CCG
	S                            P         ABE   TCT            CCT
	T                            I  Target-AID   ACA            ATA
	T                            R  Target-AID   ACA            AGA
	T                            I  Target-AID   ACC            ATC
	T                            S  Target-AID   ACC            AGC
	T                            R  Target-AID   ACG            AGG
	T                            M  Target-AID   ACG            ATG
	T                            I  Target-AID   ACT            ATT
	T                            S  Target-AID   ACT            AGT
	S                            R  Target-AID   AGC            AGG
	S                            N  Target-AID   AGC            AAC
	S                            T  Target-AID   AGC            ACC
	S                            T  Target-AID   AGT            ACT
	S                            N  Target-AID   AGT            AAT
	S                            L  Target-AID   TCA            TTA
	S                            F  Target-AID   TCC            TTC
	S                            C  Target-AID   TCC            TGC
	S                            L  Target-AID   TCG            TTG
	S                            W  Target-AID   TCG            TGG
	S                            F  Target-AID   TCT            TTT
	S                            C  Target-AID   TCT            TGT
	S can be mutated to:
	['G', 'P', 'R', 'N', 'T', 'L', 'F', 'C', 'W']
	T can be mutated to:
	['A', 'I', 'R', 'S', 'M']
	Y can be mutated to:
	['H', 'C']
	Out of total 99 sites, for 30 sites,
	guides were designed.
	Out of total 30 sites, for 13 sites,
	guides were designed for mutations to A, L or H.
	GLib-GIO-Message: Using the 'memory' GSettings backend.  Your settings will not be saved or shared with other applications.
	[2018-04-23 08:58:56,082] INFO	from pipeline.py in pipeline(..):72: Location of output data: human_test/data
	[2018-04-23 08:58:56,082] INFO	from pipeline.py in pipeline(..):73: Location of output plot: human_test/plot

## For Yeast sequences

	(01pppi)$ beditor yeast_test.json
	[2018-04-23 09:00:32,787] INFO	from pipeline.py in main(..):42: start
	[2018-04-23 09:00:34,704] INFO	from get_seq.py in din2dseq(..):124: Counts of amino acids to mutate:
	[2018-04-23 09:00:34,706] INFO	from get_seq.py in din2dseq(..):125: Y    99
	Name: aminoacid: wild-type, dtype: int64
	Possible 1 nucleotide mutations of the phospho-sites:
	           amino acid mutation method codon codon mutation
	amino acid                                                
	Y                            C    ABE   TAC            TGC
	Y                            H    ABE   TAC            CAC
	Y                            C    ABE   TAT            TGT
	Y                            H    ABE   TAT            CAT
	Y can be mutated to:
	['C', 'H']
	Out of total 99 sites, for 12 sites,
	guides were designed.
	Out of total 12 sites, for 0 sites,
	guides were designed for mutations to A, L or H.
	GLib-GIO-Message: Using the 'memory' GSettings backend.  Your settings will not be saved or shared with other applications.
	[2018-04-23 09:01:05,892] INFO	from pipeline.py in pipeline(..):72: Location of output data: yeast_test/data
	[2018-04-23 09:01:05,892] INFO	from pipeline.py in pipeline(..):73: Location of output plot: yeast_test/plot
