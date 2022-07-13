# Installation of the development version of the JWST DMS
ATOCA is not completely integrated in the official versions of the DMS pipeline software. To make the code presented here to work correctly, you need to use the `atoca_development` branch. Here are the instruction to install it:

## Option 1:
Ideally, create a new conda environment. Then do one of the following.
```
pip install git+https://github.com/AntoineDarveau/jwst.git
```

## Option 2: (allows for development)
If you are familiar with git, you can download the repository
```
git clone git@github.com:AntoineDarveau/jwst.git
cd jwst
git checkout atoca_development
pip install -e .
```
NOTE: You can also add `git@github.com:AntoineDarveau/jwst.git` as an additionnal remote repository if you already have cloned the jwst official repository.


# Example notebook to show how to use ATOCA
This [notebook](https://github.com/AntoineDarveau/atoca_demo/blob/master/Interacting_with_the_ATOCA.ipynb) is an example on how to use ATOCA within [`Extract1dStep` ](https://jwst-pipeline.readthedocs.io/en/latest/jwst/extract_1d/arguments.html) in the JWST DRS.

The inputs are described in the notebook, but in case your not able to open it, the same code is available in the file interacting_with_ATOCA.py
