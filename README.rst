=========
HFRadarPy
=========


.. .. image:: https://circleci.com/gh/rucool/HFRadarPy/tree/master.svg?style=svg
..    :target: https://circleci.com/gh/rucool/HFRadarPy/tree/master

.. .. image:: https://codecov.io/gh/rucool/hfradarpy/branch/master/graph/badge.svg
..    :target: https://codecov.io/gh/rucool/hfradarpy




Toolbox to work with High Frequency Radar (HFR) files.


* Free software: MIT license
* Documentation: https://hfradarpy.readthedocs.io.


Features
--------
* Read in CODAR, WERA and University of Hawaii HFR files
* Run quality control tests including QARTOD tests
* Calculate total vectors from radials
* Output results to CTF or NetCDF

============
Installation
============


Stable release
--------------

We recommend using Miniforge to manage your Python environments. `Download the Miniforge installer`_ for your computer architecture and operating system.
In a Terminal window, run the .sh file that was downloaded.


.. _Download the Miniforge installer: http://conda.pydata.org/miniconda.html

.. _latest release: https://conda-forge.org/download/

Add the channel, `conda-forge`_ with the following command in your Terminal:

.. code-block:: console

        conda config --add channels conda-forge

You can find out more about conda-forge from their website:

.. _conda-forge: https://conda-forge.org/

Run this command to install HFRadarPy:

.. code-block:: console

    $ mamba install hfradarpy

This method will always install the most recent stable release of HFRadarPy.

From sources
------------

The sources for HFRadarPy can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone https://github.com/rowg/HFRadarPy.git

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/rowg/hfradarpy/tarball/master

Once you have a copy of the source, you can should create a new conda/virtual environment:

Create environment
------------------

Change your current working directory to the location that you
downloaded codar_processing to.

.. code-block:: console

        $ cd ~/Downloads/hfradarpy/

Create conda environment from the included environment_dev.yml file:

.. code-block:: console

        $ conda env create -f environment_dev.yml

Once the environment is done building, you can activate the environment
by typing:

.. code-block:: console

        $ conda activate hfradarpy # OSX/Unix

Once the environment is your active environment. You can install the toolbox to that environment.

.. code-block:: console

    $ python setup.py install

You can also change directory into the root hfradarpy directory and install with the following:

.. code-block:: console

    $ pip install .

Or if you are developing new code in the toolbox, you should install this library as 'editable':

.. code-block:: console

    $ pip install --no-deps --force-reinstall --ignore-installed -e .


Running tests
-------------
After setting up your environment, you can run all of the tests, provided you install 'pytest':

.. code-block:: console

    $ pytest

.. _Github repo: https://github.com/rucool/hfradarpy
.. _tarball: https://github.com/rucool/hfradarpy/tarball/master

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
