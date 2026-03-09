# ORFfinder

A fast, parallel ORF (Open Reading Frame) finder for prokaryotic genomic sequences. Identifies protein-coding genes in FASTA sequences using the longest ORF heuristic, with support for multiple output formats.

## About ORFs

Open Reading Frames (ORFs) are continuous stretches of DNA that begin with a start codon and end with a stop codon, representing potential protein-coding genes. Key concepts:

- **Longest ORF Heuristic**: The principle that the maximum-length ORF is the most probable gene
- **Overlapping genes**: Genomic sequences where two genes overlap, potentially transcribed from opposite DNA strands
- **Translation coupling**: Adjacent genes may overlap by up to 4 nucleotides in bacteria

## Prerequisites

- Go 1.13 or later

## Installation

```bash
go build -o orffinder main.go
sudo mv orffinder /usr/local/bin/
```

Or use directly from the current directory:

```bash
./orffinder -i genome.fna -fmt fasta -o gene.fna
```

## Usage

```
orffinder -i <input.fasta> [options]
```

### Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-i` | string | (required) | Input FASTA file |
| `-min` | int | 300 | Minimum ORF length in nucleotides |
| `-fmt` | string | tsv | Output format: `tsv`, `gff3`, or `fasta` |
| `-o` | string | stdout | Output file (use `-` for stdout) |
| `-maxoverlap` | int | 4 | Max same-strand overlap (nt) between adjacent genes |
| `-translate` | flag | false | Translate ORFs to protein sequences (slower) |
| `-workers` | int | 10 | Number of parallel worker goroutines |
| `-v` | flag | false | Verbose output to stderr |

### Output Formats

- **tsv**: Tab-separated values with columns: SeqID, Strand, Frame, Start, End, Length, [Protein]
- **gff3**: GFF3 format for genome annotation tools
- **fasta**: FASTA format with nucleotide or protein sequences

### Examples

Find ORFs with default parameters and verbose output:

```bash
orffinder -i genome.fna -min 300 -fmt tsv -v
```

Generate GFF3 output with lower minimum length:

```bash
orffinder -i genome.fna -min 150 -fmt gff3 -o orfs.gff
```

Translate ORFs to protein sequences:

```bash
orffinder -i genome.fna -translate -fmt fasta -o proteins.fa
```

Use all available CPU cores for faster processing:

```bash
orffinder -i genome.fna -workers 16 -v
```
