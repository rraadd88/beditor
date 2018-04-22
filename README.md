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