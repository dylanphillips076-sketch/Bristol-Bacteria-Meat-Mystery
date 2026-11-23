#!/bin/bash
# Copy FASTA files from each outlet folder into one directory,
# adding the outlet name to each sequence header.

mkdir -p combined

for dir in data/*/; do
  outlet=$(basename "$dir")
  for file in "$dir"*.fasta; do
    name=$(basename "$file")
    sed "s/^>/>${outlet}_/" "$file" > "combined/${outlet}_${name}"
  done
done

cat combined/*.fasta > all_sequences.fasta
