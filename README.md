# Sweet orange TE insertion locus (IL) scanning and novel IL detection

## Conda environment configuration
### Download the repository
git clone https://github.com/TheLuoFengLab/Transposon-insertion-locus-simulation-and-scanning.git
mv Transposon-insertion-locus-simulation-and-scanning SWOIL  # Rename the dir to any name you want

### Install tools in TESCAN.yml
conda env create -f SWOIL/te_detection.yml

## Running IL_SCAN.sh on pair-end NGS of citrus (working for sweet orange, pummelo, and mandarins)
SWOIL/IL_SCAN.sh 


