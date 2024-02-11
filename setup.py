#!/usr/bin/env python
"""Development.

Install using 

python setup.py install

### release new version

    git commit -am "version bump";git push origin master
    python setup.py --version
    git tag -a v$(python setup.py --version) -m "upgrade";git push --tags

"""

import sys

if sys.version_info.major != 3:
    raise RuntimeError(
        f"Python 3 required. found {sys.version_info.major}"
    )
if sys.version_info.minor > 9:
    raise RuntimeError(
        f"Python <=3.9 required, because of pysam compatibility issues. found {sys.version_info.minor}"
    )

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = {
    "base": [
        # 'roux'
        "roux @ git+https://github.com/rraadd88/roux.git@master",
        ## requirements from roux
        # 'seaborn',
        # 'numpy>=1.17.3',
        # 'pandas>=0.25.3',
        # 'pyyaml>=5.1',
        # 'matplotlib>=2.2',
        "papermill",
        "biopython",
        "regex",  # expanded regex package
        "scipy",
        "pyensembl==2.2.9",
        "pysam==0.22.0",  # to convert sam to bam
        ## workflow
        "papermill",
        "argh==0.29.*",  # for command line
        # aes
        "jinja2",  # for pandas dataframe formatting
        "ipywidgets",  # for tqdm
        # 'app':[
        # "importlib-metadata==4.13.0", # to avoid celery error
        # ],
        #     ],
        # 'viz':[
        "igv-notebook==0.5.2",
        "pyfaidx==0.8.0",  ## for indexing fasta file
    ],
    "gui": [
        "mercury==2.3.7",
        "dna_features_viewer==3.1.3",
        "chrov @ git+https://github.com/rraadd88/chrov.git@main",
    ],
    ## development and maintenance
    "dev": [
        "pytest",
        "jupyter",
        "ipywidgets",
        "ipykernel",
        "black",
        "coveralls == 3.*",
        "flake8",
        "isort",
        "pytest-cov == 2.*",
        "testbook",
        # 'papermill',
        "lazydocs",  #'regex', ## docs
        # 'sphinx','recommonmark',
    ],
}
extras_require = {k: l for k, l in requirements.items() if not k == "base"}
## all: extra except dev
extras_require["all"] = [l for k, l in extras_require.items() if not k == "dev"]
### flatten
extras_require["all"] = [s for l in extras_require["all"] for s in l]
### unique
extras_require["all"] = list(set(extras_require["all"]))

## install dependency added as a submodule
# git submodule add --depth 1 --name bwa https://github.com/lh3/bwa.git beditor/bwa

import subprocess
from os.path import dirname

from setuptools.command.install import install


class InstallCommand(install):
    def run(self):
        cwd = dirname(__file__)
        com = f"cd {cwd if cwd!='' else '.'}/beditor/bwa;make"
        print(com)
        from glob import glob

        print(glob(f"{cwd}/beditor/bwa/*"))
        r = subprocess.run(com, shell=True)
        assert r.returncode == 0, r.stderr
        super().run()


# main setup
setuptools.setup(
    name="beditor",
    author="Rohan Dandage",
    version="2.0.0",
    url="https://github.com/rraadd88/beditor",
    download_url="https://github.com/rraadd88/beditor/archive/master.zip",
    description="A computational workflow for designing libraries of guide RNAs for CRISPR base editing",
    long_description="https://github.com/rraadd88/beditor/blob/master/README.md",
    keywords=["CRISPR", "genome", "biology"],
    license="General Public License v. 3",
    install_requires=requirements["base"],
    extras_require=extras_require,
    platforms="Tested on Ubuntu 22.04 64bit",
    packages=setuptools.find_packages(
        ".",
        exclude=["test", "tests", "unit", "deps", "data", "examples"],
        # include=['beditor/bwa','beditor.bwa.bwakit','beditor.data'], # ignore warning in building because otherwise there is a module not found error
    ),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "roux = beditor.run:parser.dispatch",
        ],
    },
    cmdclass={"install": InstallCommand},
    package_data={
        "data": ["beditor/data/*"],
        "bwa": ["beditor/bwa/*"],
    },
)
