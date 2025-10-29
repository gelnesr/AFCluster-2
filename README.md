# AFCluster-2

A modular reimplementation of [AF-Cluster](https://github.com/HWaymentSteele/AF_Cluster) for plug-n-play functionality, incorporation into computational workflows, and HPC/HTC deployment through ColabFold.

## Features

- **Modular Design**: Clean separation of MSA generation, clustering, and structure prediction
- **HPC-Ready**: SLURM integration for high-performance computing clusters
- **HTC-Ready**: Coming soon!
- **Modular MSA Clustering**: Additional clustering methods can be easily swapped in.

## Installation

### Prerequisites

- Python 3.8+
- ColabFold (for structure prediction)
- MMseqs2 (optional, for local MSA generation)

### Setup

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd AFCluster-2

2. **Set up ColabFold locally:**
   
   If you do not already have ColabFold isntalled, please follow the installation instructions from [localcolabfold](https://github.com/YoshitakaMo/localcolabfold):
   
   ```bash
   # For Linux
   wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
   bash install_colabbatch_linux.sh
   ```
   After setting up ColabFold, set the relative path in the `configs/afcluster.yml` file. We recommend also setting up a cache directory.

3. **For HPC deployment, run the setup script:**

    This will automatically set the path for ColabFold to your `$SCRATCH/tools` folder on HPC. Edit this line to the appropriate directory.. 
   ```bash
   bash scripts/setup_slurm.sh
   ```

## Usage

### Basic Usage

Run the pipeline on a FASTA file:

```bash
python run.py --input sequences.fasta
```

### Configuration

Modify `configs/afcluster.yml` to adjust clustering parameters:

```yaml
keyword: "MAIN"           
gap_cutoff: 0.25
random_seed: 42

cluster_method: "dbscan"

dbscan:
  min_samples: 10
  eps_val: null           
  min_eps: 3
  max_eps: 20.0
  eps_step: 0.5 

path_vars:
  PATH: "/path/to/colabfold/bin:$PATH"
  XDG_CACHE_HOME: "/path/to/cache"
  MPLCONFIGDIR: "/path/to/cache"
```

## Workflow

1. **Input Processing**: Load FASTA sequences
2. **MSA Generation**: Create multiple sequence alignments using MMseqs2
3. **Sequence Filtering**: Remove sequences with high gap content
4. **Clustering**: Group similar sequences using DBSCAN
5. **Structure Prediction**: Run ColabFold on each cluster
6. **Output**: Generate cluster-specific A3M files and predicted structures

## Output Structure

```
output/
├── sequence_id/
│   ├── sequence_id.a3m          # Generated MSA
│   ├── clusters/
│   │   ├── sequence_id_000.a3m  # Cluster 0
│   │   ├── sequence_id_001.a3m  # Cluster 1
│   │   └── ...
│   └── preds/
│       ├── sequence_id_000/
│       │   ├── s0/              # Structure prediction seed 0
│       │   ├── s1/              # Structure prediction seed 1
│       │   └── ...
│       └── ...
```

## Parameters

- `--input`: Input FASTA file (required)
- `--msa`: Pre-computed MSA file (optional)
- `--keyword`: Job identifier (default: 'MAIN')

## Dependencies

- **numpy**: Numerical computations
- **pandas**: Data manipulation
- **biopython**: Sequence processing
- **scikit-learn**: Clustering algorithms
- **tqdm**: Progress bars
- **requests**: API communication
- **pyyaml**: Configuration parsing

## HPC Requirements

- SLURM job scheduler
- GPU access (for ColabFold)
- Sufficient memory (40GB recommended)
- ColabFold installation

## License

This project is based on the original [AF-Cluster](https://github.com/HWaymentSteele/AF_Cluster) implementation. 

## Citation

If you use AFCluster-2 in your research, please cite the following works and acknowledge this implementation:

```bibtex
@article{Wayment-Steele_Ojoawo_Otten_Apitz_Pitsawong_Hömberger_Ovchinnikov_Colwell_Kern_2023,
  title={Predicting multiple conformations via sequence clustering and AlphaFold2},
  DOI={10.1038/s41586-023-06832-9},
  journal={Nature},
  author={Wayment-Steele, Hannah K. and Ojoawo, Adedolapo and Otten, Renee and Apitz, Julia M. and Pitsawong, Warintra and Hömberger, Marc and Ovchinnikov, Sergey and Colwell, Lucy and Kern, Dorothee},
  year={2023},
}
```

### Dependencies

This software builds upon the following tools and methods:

**AF-Cluster** - Wayment-Steele, H.K., Ojoawo, A., Otten, R. et al. Predicting multiple conformations via sequence clustering and AlphaFold2. *Nature* 625, 832–839 (2024). https://doi.org/10.1038/s41586-023-06832-9  
[[GitHub](https://github.com/HWaymentSteele/AF_Cluster)]

**ColabFold** - Mirdita, M., Schütze, K., Moriwaki, Y. et al. ColabFold: making protein folding accessible to all. *Nat Methods* 19, 679–682 (2022). https://doi.org/10.1038/s41592-022-01488-1  
[[GitHub](https://github.com/sokrypton/ColabFold)] | [[Local Installation](https://github.com/YoshitakaMo/localcolabfold)]

**MMseqs2** - Steinegger, M. & Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat Biotechnol* 35, 1026–1028 (2017). https://doi.org/10.1038/nbt.3981  
[[GitHub](https://github.com/soedinglab/MMseqs2)]

**AlphaFold** - Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2
