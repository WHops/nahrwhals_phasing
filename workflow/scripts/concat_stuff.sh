#!/bin/bash
set -euo pipefail
export TMPDIR=$(pwd)

# Check the number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <header.txt> <FileA.bcf> <FileB.vcf.gz> <output.vcf.gz>"
    echo "header.txt should contain the VCF header to be applied to the merged file. Use provided template as reference."
    exit 1
fi

# Assign arguments to variables for easier reference
HEADER="$1"
FILE_A="$2"
FILE_B="$3"
OUTPUT="$4"

# Check if necessary files exist
for file in "$HEADER" "$FILE_A" "$FILE_B"; do
    if [ ! -f "$file" ]; then
        echo "Error: File $file does not exist."
        exit 1
    fi
done

# Extract sample names
echo "Extracting sample names from both input files..."
bcftools query -l "$FILE_A" > samples_A.txt
bcftools query -l "$FILE_B" > samples_B.txt

# Generate sample name mapping
echo "Generating sample name mapping..."
awk 'BEGIN {
        while (getline < "samples_A.txt") {
            original = $1;
            gsub(/^GM/, "NA", $1);
            print original "\t" $1;
        }
    }
' > sample_name_mapping.txt

# Rename samples in File A
echo "Renaming samples in $FILE_A based on mapping..."
bcftools reheader -s sample_name_mapping.txt "$FILE_A" > temp_renamed.bcf

# Identify common samples
echo "Identifying common samples between the two files..."
bcftools query -l temp_renamed.bcf > samples_A_renamed.txt
awk 'NR==FNR{samples[$1]=1; next} ($1 in samples)' samples_A_renamed.txt samples_B.txt > matching_samples.txt

# Isolate relevant samples
echo "Isolating relevant samples..."
bcftools view -S matching_samples.txt temp_renamed.bcf -O b > temp_A_restricted.bcf
bcftools view -S matching_samples.txt "$FILE_B" -O b > temp_B_restricted.bcf

# Index the BCF files
echo "Indexing..."
tabix -p bcf temp_A_restricted.bcf
sleep 1
tabix -p bcf temp_B_restricted.bcf

# Concatenate the BCF files
echo "Concatenating files..."
bcftools concat temp_A_restricted.bcf temp_B_restricted.bcf -O z -a > temp_unsorted.vcf.gz

# Apply new header and sort the result
echo "Applying the new header and sorting the merged VCF..."
cat "$HEADER" <(grep -v '^##' <(zcat temp_unsorted.vcf.gz)) | bcftools sort --temp-dir $(pwd) -O z -o "$OUTPUT" -

# Cleanup intermediate files
echo "Cleaning up intermediate files..."
rm -f bcftools-sort.*
rm samples_A.txt samples_B.txt sample_name_mapping.txt temp_renamed.bcf samples_A_renamed.txt matching_samples.txt temp_A_restricted.bcf temp_B_restricted.bcf temp_unsorted.vcf.gz

echo "All done! Merged file is saved as $OUTPUT."