#     import subprocess
from beditor.lib.io_sys import runbashcmd
from os.path import exists

def test_species(host='saccharomyces_cerevisiae'):
    # clone test_beditor 
    if not exists('test_beditor'):
        runbashcmd('git clone https://github.com/rraadd88/test_beditor.git',test=True)
    else:
        runbashcmd('cd test_beditor;git pull',test=True)
    com=f'source activate beditor;cd test_beditor;python test_datasets.py'
    runbashcmd(com,test=True)
    
#     cfgp='mutation_format_nucleotide_mutation_mutations_for.yml'
#     com=f"source activate beditor;cd test_beditor/dataset_saccharomyces_cerevisiae;beditor --cfg {cfgp}"
#     runbashcmd(com,test=True)

    
test_species()
