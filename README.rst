Support Migrated
================

This repository has been migranted to `PICASO`. Please install `PICASO` and proceed to this tutorial here: 	https://natashabatalha.github.io/picaso/notebooks/workshops/SaganSchool2020/JWSTProposalTutorial.html


A Proposal Tutorial for Exoplanet Atmospheres (support migrated)
================================================================

Created for the Sagan Summer School. Huge thank you to Dr. Hannah Wakeford (Bristol) and Dr. Ehsan Gharib-Nezhad (NASA Ames) for comments and fixes. 

Setup Up 
--------

Clone the repository 

.. code-block:: bash 

	git clone https://github.com/natashabatalha/propexo.git
	cd propexo

Create conda environment from `yml` file. 

.. code-block:: bash 

	conda env create -f propexo.yml

Activate environment and open up jupyter notebook. 

.. code-block:: bash 

	conda activate propexo
	jupyter notebook

Table of Contents
-----------------

1.  What accuracy in planet properties do I need?

1.1.  Using Exoplanet Archive API to Query Confirmed Targets

2.  What tools do I need? How do I use them?

2.1.  Use PandExo to run initial constant  (𝑅𝑝/𝑅∗)2  to determine approx precision

2.2.  Use CHIMERA to determine first guess atmospheric transmission signal

2.3.  Check Exo.MAST for available data so we can validate our assumptions

2.4.  Use PICASO to determine first guess atmospheric emission signal

3.  How can I "prove" observability?

3.1  Can an atmosphere be detected: Addressing cloud concerns and quantifying statistical significance in transmission

3.2  Can an atmosphere be detected: Addressing unknown climate and quantifying statistical significance in emission

3.3  Can a specific molecule be detected?

3.4  Can any physical parameters be constrained? Information content theory for initial constraint estimates

