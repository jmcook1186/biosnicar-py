
*******
Testing
*******

BioSNICAR tests
---------------

The tests folder contains automated tests that can be run locally by navigating to the top level project directory and executing

.. code-block:: python 
   pytest ./tests


The console will show the total number of tests collected. As each test is executed its result is shown in the console as a green dot (pass) or red cross (fail). After about 5 minutes the testing will finish. 100% of tests should pass.

The tests in ./tests are also run automatically as a Github action whenever a pull request is raised against the main branch of the BioSNICAR repository. This is to help contributors and the repo owners to determine whether a PR includes breaking changes. 

The tests are all defined in `test_snicar.py`. Fixtures (features that are common to all tests) are defined in `conftest.py`. To add, remove or refine the automated tests, update these files. It is strongly recommended that these files are not changed and instead any changed are raised as issues for discussion with the repo owners because of their importance for PR quality control.

There are many tests included in `./tests` but they can be divided into two general categories: benchmarking and fuzzing.


Benchmarking
------------

The benchmark tests aim to ensure that the Python implementation in this repository is consistent with previously published versions developed and maintained by other research groups. We use the Matlab implementation of SNICAR_v3 as our benchmark model for the Toon matrix solver for spherical snow grains and the Matlab implementation of SNICAR_AD_v4 as our benchmarl for the adding-doubling solver for glacier ice. Both of these benchmark scripts were developed by Mark Flanner and Chloe Whicker at University of Michigan. We aim for wide coverage across variables to ensure consistent spectral albedo between our model and theirs.

The benchmark data was generated on the Umich linux server using the published versions of their model code. The resulting csv was saved to our BioSNICAR repository. The Python model was then run with identical config to produce another csv. These files are then compared column-wise and tested for equality. The tolerance for error between the Matlab and Python spectral albedo is 1e-8 by default but can be adjusted by updating `tolerance` in `conftest`.

The script that generates the v3 and v4 equivalent Python data is `benchmarking_driver.py`. It is **strongly recommended** to run this file after any changes to the BioSNICAR source code to ensure the benchmarking tests are using predictions from the latest version and can catch any breaking changes. This takes about 20 minutes to run and will generate new Python benchmark files. However, for rapid development the new benchmark data generation can be toggled off in conftest.

.. code-block:: python
   python ./tests/benchmarking_driver.py



Fuzzing
-------

Fuzzing aims to run the model under a wide range of config variables to check that the model runs as expected. Fuzzing checks that no combinations of variable values break the model, and attempts to catch any possible edge cases that could result from subtle bugs. The variable ranges to fuzz can be adjusted by changing the fixtures wrapping the fuzzing functions in `test_snicar.py`


Github Workflow
---------------

Any pull-requests to the main branch of BioSNICAR automatically run both the benchmarking and fuzzing tests described above in a virtual machine instance of the PR branch. Pull requests that fail any of these tests will not be merged. The complete run takes about fifteen minutes to execute - once complete the pull request will be marked with "all tests pass" so reviewers know there are no breaking changes.