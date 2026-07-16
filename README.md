[![DOI](https://zenodo.org/badge/1189627157.svg)](https://doi.org/10.5281/zenodo.21399383)

# Microbiome Analysis of Iron/Zinc-Biofortified Pearl Millet Trial

This repository contains the analysis pipeline for a randomized controlled feeding trial investigating the effects of iron/zinc-biofortified pearl millet (FeZnPM) on the gut microbiome of children in Mumbai, India (Clinical Trials.gov ID: [NCT02233764](https://clinicaltrials.gov/study/NCT02233764?term=NCT02233764&rank=1)).

## Repository Structure

```
├── analysis/     # Analysis Quarto notebooks organized by topic
├── _assets/      # Style files and images
├── manuscript/   # Manuscript files and references
├── scripts/      # Analysis and processing scripts
└── utils/        # Utility functions and helpers
```

## Analysis overview

Analysis was conducted in R 4.4.2 and package requirements and versions are available in the renv.lock file.

- Quality control and preprocessing of metagenomic data
- Classification by alignment to WoL2 reference phylogeny
- Diversity analyses (alpha and beta)
- Differential abundance testing
- Functional analysis (KEGG pathways and GO terms)
- Iron status subgroup analysis

### Docker Image

A Dockerfile has also been provided to replicate the R environment. To build the docker image locally run the following:

```bash
docker build ./ --progress=plain ./ --tag pm-microbiome &> build.log
```

```bash
dir=$(pwd)
IMAGE=biohpc_nc564/pm-microbiome-analysis

docker run -v "$dir"/:/home/workdir --rm $IMAGE quarto render workdir --log .logs/quarto-render.log --log-level info
```

## Citation

Huey, S.L., Cole, N.L., Pagani, I., González, A., Finkelstein, J.L., Haas, J.D., Udipi, S.A., Ghugre, P., Potdar, R.D., Knight, R., Mehta, S., (in preperation). Effect of a complementary feeding intervention on the gut microbiota in 12–18-month-old children.
