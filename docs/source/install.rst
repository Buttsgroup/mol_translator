Install
=======

.. _install

prerequisite
------------

As pip doesn't support the direct installation of openbabel (requires compilation), this would need to be done through conda. Recommended setup is to use miniconda.

First follow the instructions `here <https://docs.conda.io/en/latest/miniconda.html>`_ to download the corresponding version to your OS, then type the following to
create the base environment.

.. code-block:: console

    (base) $ conda create -n mol_translator python=3.9 openbabel -c conda-forge
    (base) $ conda activate mol_translator


pip
---

Simply use pip to install mol_translator once the prerequisite steps are complete

.. code-block:: console
    
    (mol_translator) pip install mol_translator

conda
-----

Installation via conda currently isn't available but it's in the works

github
------

Installation via github is also possible but this requires the most manual setup.