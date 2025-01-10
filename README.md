# Opsinth: Opsin cluster analysis with long reads

## Installlation

```bash
pip install -r requirements.txt
```

## Dependencies

Python dependencies are in `requirements.txt`.

Racon is required for polishing and can be installed from [here](https://github.com/isovic/racon). By default, the racon installation available in the PATH is used, but this can be overridden by providing a path to the `--racon` argument in the command line.

## Usage

```bash
python -m opsinth --bam <bam_file> --out <output_prefix>
```

