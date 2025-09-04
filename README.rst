==========
**Intro**
==========
- Pathannotator annotates proteins with KEGG, Flybase and Reactome 'DME' pathways. It does this through the use of `KofamScan <https://github.com/takaram/kofam_scan>`_, `KEGG API <https://www.kegg.jp/kegg/rest/keggapi.html>`_, `OrthoFinder <https://github.com/davidemms/OrthoFinder?>`_ `Flybase <https://flybase.org/>`_ and `Reactome <https://reactome.org/>`_.
- KofamScan is a gene functional annotation tool based on KEGG Orthology and hidden Markov model (HMM). It is provided by the KEGG (Kyoto Encyclopedia of Genes and Genomes) project. The online version is available here: https://www.genome.jp/tools/kofamkoala/ .
- This pipeline pulls annotation directly from the KEGG API when possible. When that isn't possible the pipeline impliments Kofamscan to identify homologous KEGG objects (KO). The pathways annotated to these KEGG objects are then transferred to the corresponding proteins in your species of interest.
- If specified, the pipeline will also provide annotations to Flybase pathways and Reactome 'DME' (*Drosophila melanogaster*). To do this the pipeline uses `OrthoFinder <https://github.com/davidemms/OrthoFinder?>`_ to identify homologous *Drosophila melanogaster* proteins for your input proteins. Flybase `metabolic pathways <http://ftp.flybase.org/releases/FB2024_06/precomputed_files/genes/metabolic_pathway_group_data_fb_2024_06.tsv.gz>`_ and `signaling pathway annotations <http://ftp.flybase.org/releases/FB2024_06/precomputed_files/genes/signaling_pathway_group_data_fb_2024_06.tsv.gz>`_ and `Reactome DME pathways <https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt>`_ are then transferred to your input proteins from these homologs.


**Where to Find Pathannotator**
=========================================

Pathannotator is provided as a Docker container for use on the command line.


- `Docker Hub <https://hub.docker.com/r/agbase/pathannotator>`_

The Dockerfile and scripts are available on `GitHub <https://github.com/AgBase/pathannotator>`_


**Getting the KOfam Databases**
===============================

The KOfam 'profiles' and 'ko_list' databases are required to run the pipeline. If you don't already have these databases the pipeline will pull them during the first run.
If you want to download them beforehand they are available from the KEGG website.

    Using wget:

    .. code-block:: bash

        wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz

        wget https://www.genome.jp/ftp/db/kofam/ko_list.gz

**TIP:**
KEGG updates their annotations approximately once a month. If your profiles and ko_list files are older than this then your annotations will not be up-to-date. Just delete them and the new versions will be downloaded with your first annotation run.

If your species of interest has been annotated by the KEGG project you can provide this tool with the corresponding KEGG species code to pull those annotations directly. If your species of interest is not listed you should choose a closely related species and use that code.
KEGG species codes can be found here: https://www.genome.jp/brite/br08611


**Help and Usage Statement**
============================
On the command line the following help statement can be displayed with 'help'.

.. code-block:: bash

    Help and Usage:
        -h to see help and usage statement
        -k KEGG species code (NA or related species code if species not in KEGG; 'help' to see this help and usage statement)
           KEGG species codes can be found here: https://www.genome.jp/brite/br08611
        -i input file (protein FASTA without header lines)
        -d (optional: default is '.') output directory (must be an existing directory; the file path should be  relative to, and inside of, your working directory)
        -f (optional: default is 'NA') 'FB' for Flybase and DME Reactome annotations, 'NA' for none
        -o outbase (file basename to use for output files)

KofamScan is used under an MIT License:

Copyright (c) 2019 Takuya Aramaki

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Flybase annotation is carried out using `OrthoFinder <https://github.com/davidemms/OrthoFinder?>`_ .

======================================
**Benchmarking**
======================================
The amount of time it takes to run this tool will vary greatly depending on several factors:
1. the number of sequences in your protein FASTA input file
2. whether your species has a KEGG species code
3. whether your input IDs are NCBI RefSeq protein IDs
4. whether you request Flybase and Reactome annotation in addition to KEGG
5. how many CPUs you have available to run the analysis

These examples below were run on the SciNet Atlas HPC system using Apptainer.

**Proteins have NCBI RefSeq IDs and that species has a KEGG code**

   +-------------------------------+------------------------+--------------------------+-------------------------------+--------------------+---------------------------+
   | **Number of input sequences** | **KEGG species code**  | **NCBI RefSeq protein**  |**Include Flybase annotation** | **Number of CPUs** | **Time to run (minutes)** |
   +-------------------------------+------------------------+--------------------------+-------------------------------+--------------------+---------------------------+
   |   22,272                      | same species           |  yes                     | yes                           | 2                  |       428                 |
   +-------------------------------+------------------------+--------------------------+-------------------------------+--------------------+---------------------------+
   |   22,272                      | same species           |  yes                     | no                            | 2                  |         < 1               |
   +-------------------------------+------------------------+--------------------------+-------------------------------+--------------------+---------------------------+
   |   22,272                      | same species           |  yes                     | yes                           | 12                 |         96                |
   +-------------------------------+------------------------+--------------------------+-------------------------------+--------------------+---------------------------+
   |   22,272                      | same species           |  yes                     | no                            | 12                 |        < 1                |
   +-------------------------------+------------------------+--------------------------+-------------------------------+--------------------+---------------------------+
   |   22,272                      | same species           |  yes                     | yes                           | 48                 |         62                |
   +-------------------------------+------------------------+--------------------------+-------------------------------+--------------------+---------------------------+
   |   22,272                      | same species           |  yes                     | no                            | 48                 |          < 1              |
   +-------------------------------+------------------------+--------------------------+-------------------------------+--------------------+---------------------------+

------------------

**Proteins have NCBI RefSeq IDs but that species doesn't have a KEGG code--using the KEGG code for a related species**

   +------------------------------+----------------------+-----------------------+------------------------------+------------------+-------------------------+
   |**Number of input sequences** | **KEGG species code**|**NCBI RefSeq protein**|**Include Flybase annotation**|**Number of CPUs**|**Time to run (minutes)**|
   +------------------------------+----------------------+-----------------------+------------------------------+------------------+-------------------------+
   |20,571                        | related species      | yes                   | yes                          | 2                |        1047             |
   +------------------------------+----------------------+-----------------------+------------------------------+------------------+-------------------------+
   |20,571                        | related species      | yes                   | no                           | 2                |         533             |
   +------------------------------+----------------------+-----------------------+------------------------------+------------------+-------------------------+
   |20,571                        | related species      | yes                   | yes                          | 12               |        192              |
   +------------------------------+----------------------+-----------------------+------------------------------+------------------+-------------------------+
   |20,571                        | related species      | yes                   | no                           | 12               |           86            |
   +------------------------------+----------------------+-----------------------+------------------------------+------------------+-------------------------+
   |20,571                        | related species      | yes                   | yes                          | 48               |        95               |
   +------------------------------+----------------------+-----------------------+------------------------------+------------------+-------------------------+
   |20,571                        | related species      | yes                   | no                           | 48               |         25              |
   +------------------------------+----------------------+-----------------------+------------------------------+------------------+-------------------------+

--------------------

**Proteins do NOT have NCBI RefSeq IDs and that species doesn't have a KEGG code--using the KEGG code for a related species**

   +-------------------------------+------------------------+--------------------------+------------------------------+-------------------+-------------------------+
   |**Number of input sequences**  | **KEGG species code**  | **NCBI RefSeq protein**  |**Include Flybase annotation**|**Number of CPUs** |**Time to run (minutes)**|
   +-------------------------------+------------------------+--------------------------+------------------------------+-------------------+-------------------------+
   |  18,330                       | related species        | no                       | yes                          | 2                 |       >12hrs            |
   +-------------------------------+------------------------+--------------------------+------------------------------+-------------------+-------------------------+
   |  18,330                       | related species        | no                       | no                           | 2                 |       199               |
   +-------------------------------+------------------------+--------------------------+------------------------------+-------------------+-------------------------+
   |  18,330                       | related species        | no                       | yes                          | 12                |       142               |
   +-------------------------------+------------------------+--------------------------+------------------------------+-------------------+-------------------------+
   |  18,330                       | related species        | no                       | no                           | 12                |        32               |
   +-------------------------------+------------------------+--------------------------+------------------------------+-------------------+-------------------------+
   |  18,330                       | related species        | no                       | yes                          | 48                |        48               |
   +-------------------------------+------------------------+--------------------------+------------------------------+-------------------+-------------------------+
   |  18,330                       | related species        | no                       | no                           | 48                |        10               |
   +-------------------------------+------------------------+--------------------------+------------------------------+-------------------+-------------------------+


======================================
**Pathannotator on the Command Line**
======================================

**Container Technologies**
===========================
Pathannotator is provided as a Docker container.

A container is a standard unit of software that packages up code and all its dependencies so the application runs quickly and reliably from one computing environment to another.

There are two major containerization technologies: **Docker** and **Apptainer (Singularity)**.

Docker containers can be run with either technology.

**Running Pathannotator using Docker**
======================================
.. admonition:: About Docker

    - Docker must be installed on the computer you wish to use for your analysis.
    - To run Docker you must have ‘root’ (admin) permissions (or use sudo).
    - Docker will run all containers as ‘root’. This makes Docker incompatible with HPC systems (see Apptainer/Singularity below).
    - Docker can be run on your local computer, a server, a cloud virtual machine etc.
    - For more information on installing Docker on other systems:  `Installing Docker <https://docs.docker.com/engine/install/>`_.


**Getting the Pathannotator container**
----------------------------------------
The Pathannotator tool is available as a Docker container on Docker Hub:
`Pathannotator container <https://hub.docker.com/r/agbase/pathannotator>`_

The container can be pulled with this command:

.. code-block:: bash

    docker pull agbase/pathannotator:4.0

.. admonition:: Remember

    You must have root permissions or use sudo, like so:

    sudo docker pull agbase/pathannotator:4.0




**Getting the Help and Usage Statement**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    sudo docker run --rm agbase/pathannotator:4.0 -h


**TIP:**

    The /workdir directory is built into this container and should be used to mount your working directory.

    The /data directory is built into this container and should be used to mount the KofamScan database files.


**Example Command**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    sudo docker run \
    --rm \
    -v /path/to/your/input/files:/workdir \
    -v /path/to/kofam/databases/:/data \
    agbase/pathannotator:4.0 \
    -k tca \
    -i GCF_031307605.1_icTriCast1.1_protein.faa \
    -d out_dir \
    -f FB \
    -o test_output

**Command Explained**
""""""""""""""""""""""

**sudo docker run:** tells docker to run

**--rm:** removes the container when the analysis has finished. The image will remain for future use.

**-v /path/to/your/input/files:/workdir:** mounts the working directory on the host machine to '/workdir' inside the container

**-v /path/to/kofam/databases/:/data:** mounts the directory with the Kofam database files (or where you want them to be stored) on the host machine to '/data' inside the container

**agbase/pathannotator:4.0:** the name of the Docker image to use

.. tip::

    All the options supplied after the image name are Pathannotator options

**-k tca:** KEGG species code for Tribolium casteneum. Can be found here: https://www.genome.jp/brite/br08611 . If your species doesn't have a code choose a closely related species or 'NA'.

**-i GCF_031307605.1_icTriCast1.1_protein.faa:** input file (protein FASTA, no header lines).

**-d out_dir:** Directory where you want the pipeline outputs to go. The directory must exist before you run the pipeline. The file path should be relative to (and inside of) your working directory.

**-f FB:** FB indicates that we want to get Flybase pathways annotations in addition to KEGG annotations.

**-o test_output:** 'test_output' will be the prefix used to name all the output files

Reference `Understanding results`_.


**Running Pathannotator using Apptainer (formerly Singularity)**
================================================================
.. admonition:: About Apptainer

    - does not require ‘root’ permissions
    - runs all containers as the user that is logged into the host machine
    - HPC systems are likely to have Apptainer installed and are unlikely to object if asked to install it (no guarantees).
    - can be run on any machine where it is installed
    - more information about `installing Apptainer <https://apptainer.org/docs-legacy>`_
    - This tool was tested using Apptainer 1.3.1

.. admonition:: HPC Job Schedulers

    Although Apptainer can be installed on any computer this documentation assumes it will be run on an HPC system. The tool was tested on a Slurm system and the job submission scripts below reflect that. Submission scripts will need to be modified for use with other job scheduler systems.

**Getting the Pathannotator container**
----------------------------------------

**Example Slurm script:**

.. code-block:: bash

    #!/bin/bash
    #SBATCH --job-name=pathannot
    #SBATCH --ntasks=8
    #SBATCH --time=2:00:00
    #SBATCH --partition=ceres
    #SBATCH --account=nal_genomics

    module load apptainer

    cd /location/where/you/want/to/save/image/file

    apptainer pull docker://agbase/pathannotator:4.0


**Running Pathannotator with Data**
------------------------------------

.. tip::

    There /workdir directory is built into this container and should be used to mount your local working directory.

    There /data directory is built into this container and should be used to mount the KOfam database files.

**Example Slurm Script**
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    #!/bin/bash
    #SBATCH --job-name=pathannot
    #SBATCH --ntasks=48
    #SBATCH --time=12:00:00
    #SBATCH --nodes=1
    #SBATCH --partition=ceres
    #SBATCH --account=nal_genomics

    module load apptainer

    cd /directory/you/want/to/work/in

    singularity run \
    -B /directory/you/want/to/work/in:/workdir \
    -B /directory/with/kofam/database/files:/data \
    /path/with/image/file/pathannotator_3.0.sif \
    -k tca \
    -i GCF_031307605.1_icTriCast1.1_protein.faa \
    -d out_dir \
    -f FB \
    -o test_output



**Command Explained**
"""""""""""""""""""""

**apptainer run:** tells Apptainer to run

**-B /directory/you/want/to/work/in:/workdir:** mounts the working directory on the host machine to '/workdir' in the container

**-B /directory/with/kofam/database/files:/data:** mounts the directory with the kofam database file (or where you want them stored) on the host machine to '/data' in the container

**/path/with/image/file/pathannotator_4.0.sif:** the name of the Apptainer image to use

.. tip::

    All the options supplied after the image name are Pathannotator options

**-k tca:** KEGG species code for Tribolium casteneum. Can be found here: https://www.genome.jp/brite/br08611 . If you species doesn't have a code choose a closely related species or 'NA'.

**-i GCF_031307605.1_icTriCast1.1_protein.faa:** input file (protein FASTA, no header lines)

**-d out_dir:** Directory where you want the outputs of the pipeline to be stored. The directory must exist before you run the pipeline. The file path should be relative to (and inside of) your working directory.

**-f FB:** FB indicates that you want Flybase pathways annotations in addition to KEGG annotations

**-o test_output:** 'test_output' will be the prefix used to name all the output files

Reference `Understanding results`_.

.. _Understanding results:

**Understanding Your Results**
==============================

The output files you can expect will differ depending on the circumstances of your run. If you are using a KEGG species code you will get both KEGG reference and KEGG species pathways. Without a KEGG code (NA) you will only get KEGG reference pathway annotations. Under all circumstances you may specify whether or not you want to receive Flybase and Reactome 'DME' pathways annotations as well. Whatever your options, the pathways will all be output into a single GMT formatted file.

**Expected output files:**
--------------------------

**Specified KEGG code**
^^^^^^^^^^^^^^^^^^^^^^^

- **test_output_KEGG_ref.tsv:** These are annotations to the KEGG reference pathways. The pathway identifiers will begin with 'map'.
- **test_output_KEGG_species.tsv:** These are annotations to the species-specific KEGG pathway. The pathway identifiers will begin with the KEGG species code.
- **test_output_flybase.tsv:** If you used the 'FB' option for Flybase pathways annotations you will get this output.
- **test_output_reactome.tsv:** If you used the 'FB' option for Flybase pathways you will get this output containg Reactome 'DME' pathways annotations.
- **test_output_all_pathways.gmt:** This file contains all of the pathways annotations (KEGG ref, KEGG species, Flybase and Reactome) in GMT format.

**'NA' as KEGG code**
^^^^^^^^^^^^^^^^^^^^^

**Expected output files:**

If you did not specify a KEGG species code (used 'NA') then no species-specific annotations file will be generated.

- **test_ouptut_KEGG_ref.tsv:** These are annotations to the KEGG reference pathways. The pathway identifiers wil begin with 'map'.
- **test_output_flybase.tsv:** If you used the 'FB' option for Flybase pathways annotations you will get this output.
- **test_output_reactome.tsv:** If you used the 'FB' option for Flybase pathways you will get this output containg Reactome 'DME' pathways annotations.
- **test_output_all_pathways.gmt:** This file contains all of the pathways annotations (KEGG ref, KEGG species, Flybase and Reactome) in GMT format.

**Expected output files:**
--------------------------

**test_ouptut_KEGG_ref.tsv:**

    +----------------------+-------------+----------------------+-------------------------------------------+
    | **Input_protein_ID** | **KEGG_KO** | **KEGG_ref_pathway** | **KEGG_ref_pathway_name**                 |
    +----------------------+-------------+----------------------+-------------------------------------------+
    |   XP_015835225.1     |  K26207     |  map04024            |    cAMP signaling pathway                 |
    +----------------------+-------------+----------------------+-------------------------------------------+
    |   XP_015835225.1     |  K26207     |  map04261            |    Adrenergic signaling in cardiomyocytes |
    +----------------------+-------------+----------------------+-------------------------------------------+
    |   XP_001813251.1     |  K01540     |  map04022            |    cGMP-PKG signaling pathway             |
    +----------------------+-------------+----------------------+-------------------------------------------+

**test_output_KEGG_species.tsv:**

    +--------------------+----------------------+----------------------+----------------------------------------------------------------------+
    |**Input_protein_ID**|**KEGG_KO**           |**KEGG_tca_pathway**  |**KEGG_tca_pathway_name**                                             |
    +--------------------+----------------------+----------------------+----------------------------------------------------------------------+
    |XP_001813251.1      |K01540                |tca04820              |Cytoskeleton in muscle cells - Tribolium castaneum (red flour beetle) |
    +--------------------+----------------------+----------------------+----------------------------------------------------------------------+
    |XP_001812480.1      |K02268                |tca00190              |Oxidative phosphorylation - Tribolium castaneum (red flour beetle)    |
    +--------------------+----------------------+----------------------+----------------------------------------------------------------------+
    |XP_008195997.1      |K04676                |tca04350              |TGF-beta signaling pathway - Tribolium castaneum (red flour beetle)   |
    +--------------------+----------------------+----------------------+----------------------------------------------------------------------+

**test_output_flybase.tsv:**

    +--------------------+----------------------+----------------------+-------------------------------------------------------+
    |**Input_protein_ID**|**Flybase_protein_ID**|**Flybase_pathway_ID**|**Flybase_pathway_name**                               |
    +--------------------+----------------------+----------------------+-------------------------------------------------------+
    |NP_001034540.1      |FBpp0077451           |FBgg0001085           |BMP Signaling Pathway Core Components                  |
    +--------------------+----------------------+----------------------+-------------------------------------------------------+
    |NP_001034503.2      |FBpp0084690           |FBgg0000904           |Insulin-like Receptor Signaling Pathway Core Components|
    +--------------------+----------------------+----------------------+-------------------------------------------------------+
    |NP_001034492.1      |FBpp0078442           |FBgg0002045           |CHITIN BIOSYNTHESIS                                    |
    +--------------------+----------------------+----------------------+-------------------------------------------------------+

**test_output_reactome.tsv:**

    +--------------------+-------------------+-----------------------+-------------------------------------------------------+
    |**Input_protein_ID**|**Uniprot_ID**     |**Reactome_pathway_ID**|**Reactome_pathway_name**                              |
    +--------------------+-------------------+-----------------------+-------------------------------------------------------+
    |XP_018221664.1      |Q9W5E1             |R-DME-1234174          |Cellular response to hypoxia                           |
    +--------------------+-------------------+-----------------------+-------------------------------------------------------+
    |XP_018222787.1      |Q9W4N8             |R-DME-194315           |Signaling by Rho GTPases                               |
    +--------------------+-------------------+-----------------------+-------------------------------------------------------+
    |XP_018222841.1      |Q7KVX1             |R-DME-162582           |Signal transduction                                    |
    +--------------------+-------------------+-----------------------+-------------------------------------------------------+

**test_output_all_pathways.gmt:**
This file does not have a header. Column1 is pathway ID, column 2 is description and columns 3-N are proteins annotated to that pathway.

    +-----------------+----------------------+-----------------+----------------+----------------+------------------------+
    |FBgg0000882      |FlyBase_pathway       |XP_018223192.1   | XP_018223020.1 | XP_018220038.1 |XP_018222736.1          |
    +-----------------+----------------------+-----------------+----------------+----------------+------------------------+
    |R-DME-109606     |Reactome_pathway      |XP_0182222056.1  |XP_018221350.1  |                |                        |
    +-----------------+----------------------+-----------------+----------------+----------------+------------------------+
    |tca00130         |KEGG_species_pathway  |XP_018219344.1   |                |                |                        |
    +-----------------+----------------------+-----------------+----------------+----------------+------------------------+
    |map00254         |KEGG_reference_pathway|XP_018220255.1   |XP_018229887.1  |                |                        |
    +-----------------+----------------------+-----------------+----------------+----------------+------------------------+


`Contact us <agbase@email.arizona.edu>`_


