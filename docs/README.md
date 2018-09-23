<div class="document">

<div class="documentwrapper">

<div class="bodywrapper">

<div class="body" role="main">

<div id="beditor" class="section">

beditor[¶](#beditor "Permalink to this headline"){.headerlink}
==============================================================

A computational workflow for designing libraries of guide RNAs for
CRISPR base editing

<div id="installation" class="section">

Installation[¶](#installation "Permalink to this headline"){.headerlink}
------------------------------------------------------------------------

1.  Create a virtual environment.

<div class="highlight-text notranslate">

<div class="highlight">

    cd beditor
    conda create -f environment.yml

</div>

</div>

2.  Activate the virtual environment.

<div class="highlight-text notranslate">

<div class="highlight">

    source activate beditor

</div>

</div>

3.  Run the analysis.

<div class="highlight-text notranslate">

<div class="highlight">

    beditor --cfg configuration.yml

</div>

</div>

</div>

<div id="configuration-file-configuration-yml" class="section">

Configuration file (configuration.yml)[¶](#configuration-file-configuration-yml "Permalink to this headline"){.headerlink}
--------------------------------------------------------------------------------------------------------------------------

It contains the all the analysis specific parameters.

Template:
<https://github.com/rraadd88/test_beditor/blob/master/common/configuration.yml>

</div>

<div id="input-format" class="section">

Input format[¶](#input-format "Permalink to this headline"){.headerlink}
------------------------------------------------------------------------

<div class="highlight-text notranslate">

<div class="highlight">

    `mutation_format` opted in configuration.yml file and corresponding columns needed in input:
    `nucleotide` : `['genome coordinate','nucleotide mutation']`.
    `aminoacid`  : `['transcript: id','aminoacid: position','amino acid mutation']`.

</div>

</div>

</div>

<div id="output-format" class="section">

Output format[¶](#output-format "Permalink to this headline"){.headerlink}
--------------------------------------------------------------------------

<div class="highlight-text notranslate">

<div class="highlight">

    `mutation_format` opted in configuration.yml file and corresponding columns needed in input:

    `nucleotide` :  ['genome coordinate','nucleotide wild-type','nucleotide mutation',]
    `nucleotide` : ['transcript: id','aminoacid: wild-type','aminoacid: position','amino acid mutation','codon: wild-type','guide: id','guide+PAM sequence','beditor score','alternate alignments count','CFD score']

</div>

</div>

Format of guide: id:

<div class="highlight-text notranslate">

<div class="highlight">

    {genomic locus}|{position}|({strategy})
    where,
    strategy= {base editor};{strand};@{distance of mutation from PAM};{PAM};{codon wild-type}:{codon mutation};{amino acid wild-type}:{amino acid mutation};

</div>

</div>

A directory by the basename of configuration file (eg. directory called
‘human’ if configuration file is ‘human.yml’) would be created in the
same folder where configuration file is located. It is referred to as
‘project directory’.

Inside a project directory there would be following folders named by
corresponding steps of analysis.

<div class="highlight-text notranslate">

<div class="highlight">

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

</div>

</div>

</div>

<div id="how-to-install-new-base-editor-or-pam" class="section">

How to install new base editor or PAM[¶](#how-to-install-new-base-editor-or-pam "Permalink to this headline"){.headerlink}
--------------------------------------------------------------------------------------------------------------------------

To list supported pams and editors

<div class="highlight-text notranslate">

<div class="highlight">

    beditor --list pams
    beditor --list editors

</div>

</div>

These lists are located in beditor/data folder and can be modified.

</div>

<div id="working-with-non-ensembl-genome" class="section">

Working with non-ensembl genome[¶](#working-with-non-ensembl-genome "Permalink to this headline"){.headerlink}
--------------------------------------------------------------------------------------------------------------

<https://github.com/openvax/pyensembl#non-ensembl-data>

</div>

<div id="how-to-analyze-test-dtasets" class="section">

How to analyze test dtasets[¶](#how-to-analyze-test-dtasets "Permalink to this headline"){.headerlink}
------------------------------------------------------------------------------------------------------

<div class="highlight-text notranslate">

<div class="highlight">

    # make the input files with mock data
    git clone https://github.com/rraadd88/test_beditor.git
    source activate beditor;cd test_beditor;python make_datasets.py
    # command to analyze mock input data of a species
    source activate beditor;cd dataset_{species name};beditor --cfg mutation_format_nucleotide_mutation_mutations_for.yml

</div>

</div>

</div>

<div id="module-beditor.lib.get_seq" class="section">

[]{#api}
API[¶](#module-beditor.lib.get_seq "Permalink to this headline"){.headerlink}
-----------------------------------------------------------------------------

 `beditor.lib.get_seq.`{.descclassname}`din2dseq`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_seq.din2dseq "Permalink to this definition"){.headerlink}

:   Wrapper for converting input data (transcript ids and positions of
    mutation) to seqeunces flanking the codon. :param cfg: configuration
    dict

<!-- -->

 `beditor.lib.get_seq.`{.descclassname}`get_seq_aminoacid`{.descname}[(]{.sig-paren}*cfg*, *din*[)]{.sig-paren}[¶](#beditor.lib.get_seq.get_seq_aminoacid "Permalink to this definition"){.headerlink}

:   Fetches sequences if mutation format is amino acid :param cfg:
    configuration dict :param din: input data :returns dsequences:
    dataframe with sequences

<!-- -->

 `beditor.lib.get_seq.`{.descclassname}`get_seq_nucleotide`{.descname}[(]{.sig-paren}*cfg*, *din*[)]{.sig-paren}[¶](#beditor.lib.get_seq.get_seq_nucleotide "Permalink to this definition"){.headerlink}

:   Fetches sequences if mutation format is nucleotide :param cfg:
    configuration dict :param din: input data :returns dsequences:
    dataframe with sequences

<!-- -->

 `beditor.lib.get_seq.`{.descclassname}`t2pmapper`{.descname}[(]{.sig-paren}*t*, *coding\_sequence\_positions*[)]{.sig-paren}[¶](#beditor.lib.get_seq.t2pmapper "Permalink to this definition"){.headerlink}

:   Maps transcript id with protein id. :param t: pyensembl transcript
    object :param t: reading frames :returns
    coding\_sequence\_positions: dataframe with mapped positions

<!-- -->

 `beditor.lib.get_seq.`{.descclassname}`tboundaries2positions`{.descname}[(]{.sig-paren}*t*[)]{.sig-paren}[¶](#beditor.lib.get_seq.tboundaries2positions "Permalink to this definition"){.headerlink}

:   Fetches positions from transcript boundaries. :param t: pyensembl
    transcript object :returns coding\_sequence\_positions: reading
    frames

[]{#module-beditor.lib.get_mutations .target}

 `beditor.lib.get_mutations.`{.descclassname}`dseq2dmutagenesis`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_mutations.dseq2dmutagenesis "Permalink to this definition"){.headerlink}

:   Generates mutagenesis strategies from identities of reference and
    mutated codons (from dseq). :param cfg: configurations from yml file

<!-- -->

 `beditor.lib.get_mutations.`{.descclassname}`filterdmutagenesis`{.descname}[(]{.sig-paren}*dmutagenesis*, *cfg*[)]{.sig-paren}[¶](#beditor.lib.get_mutations.filterdmutagenesis "Permalink to this definition"){.headerlink}

:   Filters the mutagenesis strategies by multiple options provided in
    configuration file (.yml). :param dmutagenesis: mutagenesis
    strategies (pd.DataFrame) :param cfg: configurations from yml file

<!-- -->

 `beditor.lib.get_mutations.`{.descclassname}`get_codon_table`{.descname}[(]{.sig-paren}*aa*, *tax\_id=None*[)]{.sig-paren}[¶](#beditor.lib.get_mutations.get_codon_table "Permalink to this definition"){.headerlink}

:   Gets host specific codon table. Eq:
    a\*np.exp(-(x-x0)\*\*2/(2\*sigma\*\*2)) :param aa: list of amino
    acids :param host: name of host :returns: codon table (pandas
    dataframe)

<!-- -->

 `beditor.lib.get_mutations.`{.descclassname}`get_codon_usage`{.descname}[(]{.sig-paren}*cuspp*[)]{.sig-paren}[¶](#beditor.lib.get_mutations.get_codon_usage "Permalink to this definition"){.headerlink}

:   Creates codon usage table. :param cuspp: path to cusp generated file
    :returns: codon usage table (pandas dataframe)

<!-- -->

 `beditor.lib.get_mutations.`{.descclassname}`get_possible_mutagenesis`{.descname}[(]{.sig-paren}*dcodontable*, *dcodonusage*, *BEs*, *pos\_muts*, *host*[)]{.sig-paren}[¶](#beditor.lib.get_mutations.get_possible_mutagenesis "Permalink to this definition"){.headerlink}

:   Assesses possible mutagenesis strategies, given the set of Base
    editors and positions of mutations. :param dcodontable: Codon table
    :param dcodonusage: Codon usage table :param BEs: Base editors
    (dict), see global\_vars.py :param pos\_muts: positions of mutations
    :param host: host organism :returns: possible mutagenesis strategies
    as a pandas dataframe

<!-- -->

 `beditor.lib.get_mutations.`{.descclassname}`get_submap`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_mutations.get_submap "Permalink to this definition"){.headerlink}

:   Fetches mimetic substitution map that would be used to filter
    mutagenesis strategies. Also, :param cfg: configurations from yml
    file.

[]{#module-beditor.lib.make_guides .target}

 `beditor.lib.make_guides.`{.descclassname}`dinnucleotide2dsequencesproper`{.descname}[(]{.sig-paren}*dsequences*, *dmutagenesis*, *dbug=False*[)]{.sig-paren}[¶](#beditor.lib.make_guides.dinnucleotide2dsequencesproper "Permalink to this definition"){.headerlink}

:   Makes dseqeunces dataframe of nucleotide mutation format compatible
    to guide design modules :param dsequences: dsequences dataframe
    :param dmutagenesis: dmutagenesis dataframe

<!-- -->

 `beditor.lib.make_guides.`{.descclassname}`dpam2dpam_strands`{.descname}[(]{.sig-paren}*dpam*, *pams*[)]{.sig-paren}[¶](#beditor.lib.make_guides.dpam2dpam_strands "Permalink to this definition"){.headerlink}

:   Duplicates dpam dataframe to be compatible for searching PAMs on -
    strand :param dpam: dataframe with pam information :param pams: pams
    to be used for actual designing of guides.

<!-- -->

 `beditor.lib.make_guides.`{.descclassname}`dseq2dguides`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.make_guides.dseq2dguides "Permalink to this definition"){.headerlink}

:   Wrapper around make guides function. :param cfg: configuration dict.

<!-- -->

 `beditor.lib.make_guides.`{.descclassname}`get_pam_searches`{.descname}[(]{.sig-paren}*dpam*, *seq*, *pos\_codon*, *test=False*[)]{.sig-paren}[¶](#beditor.lib.make_guides.get_pam_searches "Permalink to this definition"){.headerlink}

:   Search PAM occurance :param dpam: dataframe with PAM sequences
    :param seq: target sequence :param pos\_codon: reading frame :param
    test: debug mode on :returns dpam\_searches: dataframe with
    positions of pams

<!-- -->

 `beditor.lib.make_guides.`{.descclassname}`guide2dpositions`{.descname}[(]{.sig-paren}*x*, *dbug=False*[)]{.sig-paren}[¶](#beditor.lib.make_guides.guide2dpositions "Permalink to this definition"){.headerlink}

:   Get positions of guides relative to the target site and PAM sequence
    Note: Index and flank sequence based indexing are 0-based Distances
    and positions from pam are 1-based :param x: lambda section of
    dguides dataframe

<!-- -->

 `beditor.lib.make_guides.`{.descclassname}`make_guides`{.descname}[(]{.sig-paren}*cfg*, *dseq*, *dmutagenesis*, *dpam*, *test=False*, *dbug=False*[)]{.sig-paren}[¶](#beditor.lib.make_guides.make_guides "Permalink to this definition"){.headerlink}

:   Wrapper around submodules that design guides by 1. searching all PAM
    sequences on ‘both’ the strands, 2. filtering guides by all possible
    strategies (given in dmutagenesis) e.g. activity window, Finally
    generates a table. :param cfg: configuration dict :param dseq:
    dsequences dataframe :param dmutagenesis: dmutagenesis dataframe
    :param dpam: dpam dataframe :param test: debug mode on :param dbug:
    more verbose

[]{#module-beditor.lib.get_specificity .target}

 `beditor.lib.get_specificity.`{.descclassname}`alignmentbed2dalignedfasta`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.alignmentbed2dalignedfasta "Permalink to this definition"){.headerlink}

:   Get sequences in FASTA format from BED file step\#5 :param cfg:
    configuration dict

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dalignbed2annotationsbed`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dalignbed2annotationsbed "Permalink to this definition"){.headerlink}

:   Get annotations from the aligned BED file step\#3 :param cfg:
    configuration dict

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dalignbed2dalignbedguides`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dalignbed2dalignbedguides "Permalink to this definition"){.headerlink}

:   Get guide seqeunces from the BED file step\#4 :param cfg:
    configuration dict

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dalignbed2dalignbedguidesseq`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dalignbed2dalignbedguidesseq "Permalink to this definition"){.headerlink}

:   Get sequences from BED file step\#6 :param cfg: configuration dict

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dalignbedannot2daggbyguide`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dalignbedannot2daggbyguide "Permalink to this definition"){.headerlink}

:   Aggregate annotations per alignment to annotations per guide.
    step\#10 :param cfg: configuration dict

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dalignbedguidesseq2dalignbedstats`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dalignbedguidesseq2dalignbedstats "Permalink to this definition"){.headerlink}

:   Gets scores for guides step\#7 :param cfg: configuration dict

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dannots2dalignbed2dannotsagg`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dannots2dalignbed2dannotsagg "Permalink to this definition"){.headerlink}

:   Aggregate annotations per guide step\#8 :param cfg: configuration
    dict

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dannotsagg2dannots2dalignbedannot`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dannotsagg2dannots2dalignbedannot "Permalink to this definition"){.headerlink}

:   Map aggregated annotations to guides step\#9 :param cfg:
    configuration dict

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dguides2guidessam`{.descname}[(]{.sig-paren}*cfg*, *dguides*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dguides2guidessam "Permalink to this definition"){.headerlink}

:   Aligns guides to genome and gets SAM file step\#1 :param cfg:
    configuration dict :param dguides: dataframe of guides

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`dguides2offtargets`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.dguides2offtargets "Permalink to this definition"){.headerlink}

:   All the processes in offtarget detection are here. :param cfg:
    Configuration settings provided in .yml file

<!-- -->

 `beditor.lib.get_specificity.`{.descclassname}`guidessam2dalignbed`{.descname}[(]{.sig-paren}*cfg*[)]{.sig-paren}[¶](#beditor.lib.get_specificity.guidessam2dalignbed "Permalink to this definition"){.headerlink}

:   Processes SAM file to get the genomic coordinates in BED format
    step\#2 :param cfg: configuration dict

</div>

</div>

</div>

</div>

</div>

<div class="sphinxsidebar" role="navigation"
aria-label="main navigation">

<div class="sphinxsidebarwrapper">

[beditor](#) {#beditor-1 .logo}
============

### Navigation

<div class="relations">

### Related Topics

-   [Documentation overview](#)

</div>

<div id="searchbox" style="display: none" role="search">

### Quick search

<div class="searchformwrapper">

</div>

</div>

</div>

</div>

<div class="clearer">

</div>

</div>

<div class="footer">

©2018, Rohan Dandage. | Powered by [Sphinx
1.8.0](http://sphinx-doc.org/) & [Alabaster
0.7.11](https://github.com/bitprophet/alabaster) | [Page
source](_sources/README.rst.txt)

</div>
