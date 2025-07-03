| Usage | Requirement |
| :--- | :--- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v24.10.4-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | [![Dependencies: Apptainer Version](https://img.shields.io/badge/Apptainer-v1.3.5-blue?style=plastic)](https://github.com/apptainer/apptainer) |
| | [![Dependencies: Graphviz Version](https://img.shields.io/badge/Graphviz-v2.42.3-blue?style=plastic)](https://www.graphviz.org/download/) |

<br /><br />
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [WARNING](#warning)
   - [CONTENT](#content)
   - [INPUT](#input)
   - [HOW TO RUN](#how-to-run)
   - [OUTPUT](#output)
   - [VERSIONS](#versions)
   - [LICENCE](#licence)
   - [CITATION](#citation)
   - [CREDITS](#credits)
   - [ACKNOWLEDGEMENTS](#Acknowledgements)
   - [WHAT'S NEW IN](#what's-new-in)

<br /><br />
## AIM

Comparative analysis of methylation in E. coli using Hifi PacBio Revio long reads.

<br /><br />
## WARNINGS

Can be only used for specific experimental designs performed by Loot et al.

<br /><br />
## CONTENT

| Files and folder | Description |
| :--- | :--- |
| **main.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. |
| **nextflow.config** | Parameter settings for the *main.nf* file. Users have to open this file, set the desired settings and save these modifications before execution. Of note, this configuration file is systematically saved in the reports folder (see [below](#output)) during each execution, to save the parameter settings. |
| **bin folder** | Contains files required by the *main.nf* file. |
| **Licence.txt** | Licence of the release. |


<br /><br />
## INPUT

| Required files |
| :--- |
| A zipped or not hifi_reads.bam file. |
| A reference genome fasta file. |


<br />

The dataset used in the *nextflow.config* file, as an example, is available at https://zenodo.org/records/15798125/files/W621.hifi_reads.bam.


<br /><br />
## HOW TO RUN

### 1. Prerequisite

Installation of:<br />
[nextflow DSL2](https://gael-millot.github.io/protocols/docs/Protocol%20152-rev0%20DSL2.html#_Toc159933761). Please, use the version indicated above.<br />
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu.<br />
[Apptainer](https://gael-millot.github.io/protocols/docs/Protocol%20135-rev0%20APPTAINER.html#_Toc160091693).<br />
<br />

Optional installation (to avoid reccurent message) of:<br />
[Gocryptfs](https://github.com/rfjakob/gocryptfs), `sudo apt install gocryptfs` for Linux ubuntu.<br /> 

<br />

### 2. Local running (personal computer)

#### 2.1. *main.nf* file in the personal computer

- Mount a server if required:

<pre>
DRIVE="Z" # change the letter to fit the correct drive
sudo mkdir /mnt/share
sudo mount -t drvfs $DRIVE: /mnt/share
</pre>

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like:
<pre>
Launching `main.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
</pre>

- Run the following command from where the *main.nf* and *nextflow.config* files are (example: \\wsl$\Ubuntu-20.04\home\gael):

<pre>
nextflow run main.nf -c nextflow.config
</pre>

with -c to specify the name of the config file used.

<br /><br />
#### 2.3. *main.nf* file in the public gitlab repository

Run the following command from where you want the results:

<pre>
nextflow run -hub pasteur gmillot/19583_loot -r v1.0.0
</pre>

<br /><br />
### 3. Distant running (example with the Pasteur cluster)

#### 3.1. Pre-execution

Copy-paste this after having modified the EXEC_PATH variable:

<pre>
EXEC_PATH="/pasteur/helix/projects/BioIT/gmillot/19583_loot" # where the bin folder of 19583_loot is located (by default, the same path as for the main.nf file)
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export APP_CONF=apptainer/1.3.5
export APP_CONF_AFTER=bin/apptainer # on maestro
export GIT_CONF=git/2.39.1
export GIT_CONF_AFTER=bin/git # on maestro
export GRAPHVIZ_CONF=graphviz/2.42.3
export GRAPHVIZ_CONF_AFTER=bin/graphviz # on maestro
export GRAALVM_CONF=graalvm/ce-java17-22.3.1 # required for nextflow
export GRAALVM_CONF_AFTER=bin/graalvm # on maestro
export NEXTFLOW_CONF=nextflow/24.10.3
export NEXTFLOW_CONF_AFTER=bin/nextflow # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${APP_CONF}/${APP_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER},${CONF_BEFORE}/${GRAPHVIZ_CONF}/${GRAPHVIZ_CONF_AFTER},${CONF_BEFORE}/${GRAALVM_CONF}/${GRAALVM_CONF_AFTER},${CONF_BEFORE}/${NEXTFLOW_CONF}/${NEXTFLOW_CONF_AFTER}"
cd ${EXEC_PATH}
chmod 755 ${EXEC_PATH}/bin/*.* # not required if no bin folder
module load ${JAVA_CONF} ${APP_CONF} ${GIT_CONF} ${GRAPHVIZ_CONF} ${GRAALVM_CONF} ${NEXTFLOW_CONF}
</pre>

<br /><br />
#### 3.2. *main.nf* file in a cluster folder

Modify the second line of the code below, and run from where the *main.nf* and *nextflow.config* files are (which has been set thanks to the EXEC_PATH variable above):

<pre>
HOME_INI=$HOME
HOME="${HELIXHOME}/19583_loot/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/19583_loot/, for instance. See NFX_HOME in the nextflow software script
nextflow run main.nf -c nextflow.config
HOME=$HOME_INI
</pre>

<br /><br />
#### 3.3. *main.nf* file in the public gitlab repository

Modify the first and third lines of the code below, and run (results will be where the EXEC_PATH variable has been set above):

<pre>
VERSION="v1.0"
HOME_INI=$HOME
HOME="${HELIXHOME}/19583_loot/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/19583_loot/, for instance. See NFX_HOME in the nextflow software script
nextflow run -hub pasteur gmillot/19583_loot -r $VERSION -c $HOME/nextflow.config
HOME=$HOME_INI
</pre>

<br /><br />
### 4. Error messages and solutions

#### Message 1
```
Unknown error accessing project `gmillot/19583_loot` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/19583_loot
```

Purge using:
<pre>
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
</pre>

#### Message 2
```
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2F19583_loot
```

Contact Gael Millot (distant repository is not public).

#### Message 3

```
permission denied
```

Use chmod to change the user rights. Example linked to files in the bin folder: 
```
chmod 755 bin/*.*
```


<br /><br />

## OUTPUT

By default, all the results are returned in a *result* folder where the *main.nf* executed file is located (created if does not exist). This can be changed using the *out_path_ini* parameter of the *nextflow.config* file. By default, each execution produces a new folder named *19583_loot_\<ID\>*, created inside the *result* folder and containing all the outputs of the execution. The name of the folder can be changed using the *result_folder_name* parameter of the *nextflow.config* file. The new name file will be followed by an \<ID\> in all cases.
<br /><br />
An example of results obtained with the dataset is present at this address: https://zenodo.org/records/15132203/files/19583_loot_1743690584.zip.

<br /><br />
Mandatory elements:
<br /><br />
| 19583_loot_<UNIQUE_ID> folders and files | Description |
| --- | --- |
| **report.html** | Report of the analysis. |
| **reports** | Folder containing all the reports of the different processes as well as the **nextflow.config** file used. |
| **files** | Folder containing some of the output files of the processes. |



<br /><br />
## VERSIONS


The different releases are tagged [here](https://gitlab.pasteur.fr/gmillot/19583_loot/-/tags).

<br /><br />
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchandability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses or in the Licence.txt attached file.


<br /><br />
## CITATION

Not yet published


<br /><br />
## CREDITS

[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Institut Pasteur, Université Paris Cité, Bioinformatics and Biostatistics Hub, 75015 Paris, France

Brice Raffestin, HPC Core Facility, Institut Pasteur, Université Paris Cité, Bioinformatics and Biostatistics Hub, 75015 Paris, France

<br /><br />
## ACKNOWLEDGEMENTS


The developers & maintainers of the mentioned softwares and packages, including:

- [R](https://www.r-project.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [Nextflow](https://www.nextflow.io/)
- [Apptainer](https://apptainer.org/)
- [Docker](https://www.docker.com/)
- [Gitlab](https://about.gitlab.com/)
- [Bash](https://www.gnu.org/software/bash/)
- [Ubuntu](https://ubuntu.com/)

Special acknowledgement to:

- [PacBio](https://www.pacb.com/) and its [support team](support@pacb.com) for the release of [SMRT Link](https://www.pacb.com/wp-content/uploads/SMRT-Link-software-installation-guide-v25.2.pdf).

<br /><br />
## WHAT'S NEW IN

#### v1.0

- Everything.



