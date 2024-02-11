#!usr/bin/python
"""Command-line options"""
import argh
import logging

# from os.path import exists
from pathlib import Path
from roux.workflow.task import run_tasks
from beditor.lib.utils import get_src_path

def validate_params(
    parameters: dict,
) -> bool:
    """Validate the parameters.

    Args:
        parameters (dict): parameters

    Returns:
        bool: whther the parameters are valid or not
    """
    for k, v in parameters.items():
        if (
            k.endswith("_path")
            and v is not None
            and k not in ["output_path", "output_dir_path"]
        ):
            assert Path(v).exists(), f"Error: Path not found {v} ({k})"
        if k in ["method"]:
            assert (
                parameters[k] is not None
            ), "Error: Editor/method is a requred parameter"
    assert (
        sum([parameters[k] is not None for k in ["species_name", "release"]]) == 2
        or sum(
            [
                parameters[k] is not None
                for k in ["genome_path", "gtf_path", "transcript_path", "protein_path"]
            ]
        )
        >= 2
    ), f"Either Ensembl genome information ({[parameters[k] for k in ['species_name','release']]}) or the paths to the files ({[parameters[k] for k in ['genome_path','gtf_path','transcript_path','protein_path']]}) are needed (not both). "
    return True


def cli(
    editor: str = None,
    mutations_path: str = None,
    output_dir_path: str = None,
    species: str = None,
    ensembl_release: int = None,
    genome_path: str = None,
    gtf_path: str = None,  #'../examples/inputs/ann.gtf',
    rna_path: str = None,  # "../examples/inputs/RNA.fa",
    prt_path: str = None,  # "../examples/inputs/Protein.fa",
    search_window: int = None,
    not_be: bool = False,
    # or
    config_path: str = None,
    wd_path: str = None,
    threads: int = 1,
    kernel_name: str = "beditor",
    verbose="WARNING",
    igv_path_prefix=None,
    ext: str = None,
    force: bool = False,
    dbug: bool = False,
    skip=None,
    **kws,
):
    """beditor command-line (CLI) 

    Args:
        editor (str, optional): base-editing method, available methods can be listed using command: 'beditor resources'. Defaults to None.
        mutations_path (str, optional): path to the mutation file, the format of which is available at https://github.com/rraadd88/beditor/README.md#Input-format. Defaults to None.
        output_dir_path (str, optional): path to the directory where the outputs should be saved. Defaults to None.
        species (str, optional): species name. Defaults to None.
        ensembl_release (int, optional): ensemble release number. Defaults to None.
        genome_path (str, optional): path to the genome file, which is not available on Ensembl. Defaults to None.
        gtf_path (str, optional): path to the gene annotations file, which is not available on Ensembl. Defaults to None.
        rna_path (str, optional): path to the transcript sequences file, which is not available on Ensembl. Defaults to None.
        prt_path (str, optional): path to the protein sequences file, which is not available on Ensembl. Defaults to None.
        search_window (int, optional): number of bases to search on either side of a target, if not specified, it is inferred by beditor. Defaults to None.
        not_be (bool, optional): do not process as a base editor. Defaults to False.
        config_path (str, optional): path to the configuration file. Defaults to None.
        wd_path (str, optional): path to the working directory. Defaults to None.
        threads (int, optional): number of threads. Defaults to 1.
        kernel_name (str, optional): name of the jupyter kernel. Defaults to "beditor".
        verbose (str, optional): verbose, logging levels: DEBUG > INFO > WARNING > ERROR (default) > CRITICAL. Defaults to "WARNING".
        igv_path_prefix (_type_, optional): prefix to be added to the IGV url. Defaults to None.
        ext (str, optional): file extensions of the output tables. Defaults to None.
        force (bool, optional): overwrite the outputs of they exist. Defaults to False.
        dbug (bool, optional): debug mode (developer). Defaults to False.
        skip (_type_, optional): skip sections of the workflow (developer). Defaults to None.
        
    Examples:
        beditor cli -c inputs/mutations/protein/positions.yml
        
    Notes:
        Required parameters for a run:
            editor
            mutations_path
            output_dir_path

            or

            config_path
    """
    from datetime import datetime
    from roux.lib.io import to_dict

    _start_time = datetime.now()

    if editor is None and config_path is None:
        logging.error(
            "Required parameters not provided. Use --help to get more information."
        )
        return
    if threads is None:
        threads = 1
    ## setting verbose
    logging.basicConfig(
        level=getattr(logging, verbose),
        format="[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s",
        force=True,
    )

    if config_path is not None:
        from roux.lib.io import read_dict

        parameters = read_dict(config_path)["run"]
        ## over-write with commandline arguments
        if editor is not None:
            logging.warning("using command-line argument for editor")
            parameters["method"] = editor
        if mutations_path is not None:
            logging.warning("using command-line argument for input_path")
            parameters["input_path"] = mutations_path
        if output_dir_path is not None:
            logging.warning("using command-line argument for output_dir_path")
            parameters["output_dir_path"] = output_dir_path
        if ext is not None:
            logging.warning("using command-line argument for ext")
            parameters["output_ext"] = ext
    else:
        parameters = dict(
            method=editor,  #'Cas12a-BE', # one method for one run
            input_path=mutations_path,  #'../examples/inputs/mutations.tsv',
            output_dir_path=output_dir_path,
            search_window=search_window,  ## bases left and right of the target to define the region to design guides in
            ## for species registered in pyensembl
            species_name=species,
            release=ensembl_release,
            ## for non-registered species
            genome_path=genome_path,  #'../examples/inputs/dna.fa',
            gtf_path=gtf_path,  #'../examples/inputs/ann.gtf',
            transcript_path=rna_path,  # "../examples/inputs/RNA.fa",
            protein_path=prt_path,  # "../examples/inputs/Protein.fa",
            not_be=not_be,  # True, # non-base editing applications (skips filtering by the editable base)
            # threads=threads,
            # force=force,
            # verbose=verbose=='DEBUG',
        )
    ## overwrite with command line parameters
    parameters["threads"] = threads
    parameters["force"] = force
    parameters["verbose"] = verbose == "DEBUG"
    parameters["dbug"] = dbug
    parameters["igv_path_prefix"] = igv_path_prefix
    parameters["wd_path"] = wd_path

    if "output_path" not in parameters:
        parameters[
            "output_path"
        ] = f"{parameters['output_dir_path']}/output.tsv"  #'../examples/outputs/05_small.tsv',
        del parameters["output_dir_path"]
    if parameters["output_path"] is None:
        parameters[
            "output_path"
        ] = f"{parameters['output_dir_path']}/output.tsv"  #'../examples/outputs/05_small.tsv',
        del parameters["output_dir_path"]
    
    if parameters["wd_path"] is None:
        parameters["wd_path"] = str(Path("./").absolute())

    logging.debug(f"Parameters provided: {parameters}")

    ## validate the parameters
    validate_params(parameters)
    ## validate input mutations
    ## validate method
    if skip is not None:
        if isinstance(skip, str):
            skip = skip.split(",")
        assert isinstance(skip, list)
        logging.warning(f"skip = {skip}")
    else:
        skip = []

    ## output config
    config = {"run": parameters}
    del parameters
    ## config used henceforth
    config_out_path = f"{Path(config['run']['output_path']).parent}/config.yaml"
    to_dict(config, config_out_path)

    paths = run_tasks(
        input_notebook_path=f"{get_src_path()}/workflow.ipynb",
        kernel=kernel_name,
        parameters_list=[config["run"]],
        to_filter_nbby_patterns_kws=dict(
            patterns=[
                " app",
                "Demo",  # do not remove in app
            ]
            + skip
        ),
        force=force,
        out_paths=True,
        **kws,
    )

    from beditor.lib.utils import log_time_elapsed

    print(f"Time taken is saved at       : {log_time_elapsed(_start_time)}")
    print(f"The configuration is saved at: {config_out_path}")
    if len(paths) != 0:
        config["report path"] = paths[0]
        print(f"The report is saved at       : {config['report path']}")
    print(f"The output is saved at       : {config['run']['output_path']}")


def gui():
    logging.warning("GUI is recommended for: 1. small genomes,")
    logging.warning("                        2. small number of input mutations, and ")
    logging.warning(
        "                        3. for genomes that have been already installed by beditor's command line interface (cli)."
    )

    port = 8000
    print(f"Direct link to the beditor gui: http://127.0.0.1:{port}/app/beditor")
    nb_path = f"{get_src_path()}/beditor.ipynb"
    import subprocess

    com = f"mercury run clear;mercury run {nb_path} --verbose"
    return subprocess.run(com, shell=True)

def resources():
    methods_path = f"{get_src_path()}/data/methods.tsv"
    from roux.lib.io import read_table
    dbepams = read_table(methods_path, params=dict(sep="\t", keep_default_na=False))
    print("Supported methods:")
    print(dbepams.drop(['PAM position','distance of mutation from PAM: minimum','distance of mutation from PAM: maximum'],axis=1))
    from beditor.lib.utils import get_annots_installed
    print("") 
    print("Installed genomes:") 
    print(get_annots_installed())

from beditor.lib.io import to_fasta_index,to_2bit

parser = argh.ArghParser()
parser.add_commands(
    [
        cli,
        gui,
        resources,
        to_2bit,
        to_fasta_index
    ]
)

if __name__ == "__main__":
    parser.dispatch()
