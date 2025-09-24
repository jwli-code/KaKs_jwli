# KaKs_jwli 软件安装
### conda 环境
`micromamba create -n kaks && micromamba activate kaks && micromamba install -y muscle trimal gblocks pal2nal kakscalculator2`
### Git 软件包
`git clone https://github.com/jwli-code/KaKs_jwli.git`
## Description
This tool calculates Ka/Ks ratios for gene pairs using:
- MUSCLE for sequence alignment
- PAL2NAL for protein-to-nucleotide alignment conversion
- TrimAl/Gblocks for sequence trimming
- KaKs_Calculator for Ka/Ks ratio calculation

## Required Arguments

| Option | Short Form | Description |
|--------|------------|-------------|
| `--genepair` | `-p` | Gene pair list file (tab-delimited format) |
| `--cds1` | `-1` | CDS sequences for species 1 |
| `--cds2` | `-2` | CDS sequences for species 2 |
| `--pep1` | `-a` | Protein sequences for species 1 |
| `--pep2` | `-b` | Protein sequences for species 2 |
| `--output` | `-o` | Output directory for results |

## Optional Arguments

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Show this help message and exit |

## Example Usage
