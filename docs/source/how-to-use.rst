**********
How to Use
**********

This document explains how to use BioSNICAR.

Installing Environment/Dependencies
-----------------------------------

If you do not have Python installed, download Python >3.6. It is recommended to use a fresh environment using conda or venv. Once activated, install the project dependencies with:

.. code-block::python3
    pip install -r requirements.txt

Finally, if you do not wish to install anything on your computer, but you use VSCode and Docker, then you can use the devcontainer config provided to run this code in a remote container. This requires the "remote containers" extension to be added to VSCode. Further instructions are available here: https://code.visualstudio.com/docs/remote/containers


Running the model
-----------------

The model driver and all the core source code can be found in `/src`. From the top level directory (`~/BioSNICAR_GO_PY`), run:

`python ./src/snicar_driver.py`

This will run the model with all the default settings. The user will see a list of output values printed to the console and a spectral albedo plot appear in a separate window. The code can also be run in an interactive session (Jupyter/iPython) in which case the relevant data and figure will appear in the interactive console. 

Most users will want to experiment with changing input parameters. This is achieved by adjusting the values in the config file `inputs.yaml`. The nature of each parameter is described in in-line annotations to guide the user. Invalid combinations of values will be rejected by our error-checking code. Most users should have no reason to modify any other file in this repository except for the those in `inputs.yaml`.

More complex applications of the model code, for example model inversions, field/model comparisons etc are included under `/experiments`, with details provided in that module's own README.

We have also maintained a separate version of the BioSNICAR codebase that uses a "functional" prorgamming style rather than the object-oriented approach taken here. We refer to this as BioSNICAR Classic and it is available in the `classic` branch of this repository. it might be useful for people already familiar with FORTRAN or Matblab implementations from previous literature. The two branches are entirely equivalent int heir simulations but very different in their programmign style. The object oriented approach is preferred because it is more Pythonic, more flexible and easier to debug.

Choosing Inputs
------------------
It is straightforward to adjust the model configuration by updating the values in `inputs.yaml`. However there is a lot of nuance to setting up the model to provide realistic simulations, and the meaning of the various parameters is not always obvious. for this reason we have put together a guide. Please refer to `inputs.md` for explanations of each parameter. 