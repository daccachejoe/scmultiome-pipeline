#!/bin/bash

# install the required packages to make scmultiome executable in the desired directory
git clone https://bitbucket.org/djhshih/argparser.git
cd argparser
Rscript install.R
R CMD INSTALL .

chmod +x test.R
chmod +x multiome-processing.R