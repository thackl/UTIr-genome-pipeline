#!/usr/bin/env bash

if [ $# -ne 2 ]; then
    echo "Create a relative link for a given file from any working directory"
    echo "Usage: $0 path/to/file.txt path/to/dir/[link]";
    exit;
fi;

input="$1"
output="$2"

# Determine if output is a directory
if [ -d "$output" ]; then
  # Output is a directory: symlink goes into directory, preserving input filename
  output_path="$output/$(basename "$input")"
  output_dir="$output"
else
  # Output is a file path: use as-is
  output_path="$output"
  output_dir=$(dirname "$output")
fi

# Calculate relative path from output directory to input file
rel_path=$(realpath --relative-to="$output_dir" "$input")

# Create symlink
ln --force -s "$rel_path" "$output_path"

echo "Linking: $output_path -> $rel_path"
