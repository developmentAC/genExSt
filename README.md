![logo](graphics/genExST_logo.png)

## GenExSt (as in, "Gene-next")

A Tool to Identify Correlation of Gene Expression after Normalization with Housekeeping Genes

[![MIT Licence](https://img.shields.io/bower/l/bootstrap)](https://opensource.org/licenses/MIT)

[![Built with Streamlit](https://img.shields.io/badge/built%20with-Streamlit-09a3d5.svg)](https://www.streamlit.io/)

[![BLM](https://img.shields.io/badge/BlackLivesMatter-yellow)](https://blacklivesmatter.com/)

Stay tuned; more documentation is coming. :-)

---

## Table of contents

* [Contact](#contact)
* [Overview](#overview)
* [WebSites](#websites)
* [Externals](#external)
* [Docker](#docker)
* [References](#references)

---

## Contact
 - Oliver Bonham-Carter,
 - Date: 29 April 2021
 - email: obonhamcarter@allegheny.edu

## WebSites
 - [Department of Computer Science](https://www.cs.allegheny.edu/sites/obonhamcarter/)
 - [Allegheny College](https://allegheny.edu/)

## Externals

Here, the reader will find links to external resources for the tool.

 - GitHub link: https://github.com/developmentAC/genExSt

 - YouTube video discussion and demonstration: https://www.youtube.com/watch?v=0p6iL9wD_PU

## Overview

A tool to automate gene expression normalization (via housekeeping genes) for
studying potential correlations across a wide and diverse set of gene expressions.

### Docker

So that the user can be sure that all necessary libraries have been conveniently
installed and are present such as Streamlit (https://www.streamlit.io/), Plotly
(https://plotly.com/) and SciPy (https://www.scipy.org/), it is recommended that
Docker Desktop be utilized for running GenExSt.

Once Docker Desktop has been installed from (https://www.docker.com/), the
Docker container can be built using the OS-specific commands which are shown below.

#### OS-specific scripts to build and run containers

The following bash scripts simplify building the container.

| OS  | Building  | Running  |
|---|---|---|
| MacOS  		|  `./build_macOS.sh` |  `./run_macOS.sh` |
| Linux   	|  `./build_linux.sh` | `./run_linux.sh`  |
| Windows 	|  `build_win.bat` 		|  `run_win.bat` |

These files may be found in the directory, `dockerRunScripts/` and the builder require a copy of `Dockerfile` to run. The `Dockerfile` is found in the main directory and so it is recommended that the user stay in the main and enter the command, ` sh ./dockerRunScripts/build_macOS.sh` or similar to build a container. To run the container, type the command ` sh ./dockerRunScripts/run_macOS.sh`. Us an equivalent command for each OS.

Please note that you may be required to enter your password twice,
depending on your machine. The first time you enter your password
will be to build and initialize the Docker container. The second
time you enter your password will be to change ownership of your
output files from `root` to `$USER` once you exit the container.
These commands are maintained in the run script file.

## References

This tool has been published at the FICC2021 conference in the below article.

Oliver Bonham-Carter and Yee Mon Thu, "GenExSt: A Tool to Identify Correlation
of Gene Expression after Normalization with Housekeeping Genes.", Future
of Information and Communication Conference, Advances in Information and
Communication, Proceedings of the 2021 Future of Information and Communication
Conference (FICC), Volume 2.

Note: Once the official reference is known, it will be updated here.

---

## A Work In Progress

GenExSt is a work-in-progress. Tests for the code will come soon and I will
continue to update the repository with updates. If you would like to contribute
to this project, __then please do!__ For instance, if you see some low-hanging
fruit or task that you could easily complete, that could add value to the project,
then I would love to have your insight.

Otherwise, please create an Issue for bugs or errors. Since I am a teaching
faculty member at Allegheny College, I may not have all the time necessary
to quickly fix the bugs and so I would be very happy to have any help that
I can get from the OpenSource community for any technological insight. Much
thanks in advance. I hope that this project helps you find the knowledge from
PubMed that you require for your project. :-)
