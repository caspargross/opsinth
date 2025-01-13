# Opsinth: Opsin cluster analysis with long reads

## Installlation

clone this package from GitHub
```
git clone 
```

### Option A

Install python dependencies with pipx or with pip into a new virtualenv
```
pipx -f requirements.txt
```

We also need racon for polishing. Install it from source [here](https://github.com/isovic/racon) and add the executable to your PATH environment. The path to the racon executable can also be provided with the `--racon-path` command line argument. 

### Option B: Install with conda

Alternatively you can also install racon together with all required python dependencies via conda
```
conda env create -f conda_env.yml -n opsinth
conda activate opsinth
```
## Usage

```bash
python -m opsinth --bam <bam_file> --out <output_prefix>
```

