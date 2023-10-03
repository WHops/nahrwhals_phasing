#!/bin/bash

# Check the number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <FileA.bcf> <FileB.vcf.gz> <output.vcf.gz>"
    exit 1
fi

# Assign arguments to variables for easier reference
FILE_A="$1"
FILE_B="$2"
OUTPUT="$3"

# Check if input files have index
if [ ! -f "${FILE_A}.csi" ] && [ ! -f "${FILE_A}.csi" ]; then
    echo "Error: Index not found for $FILE_A"
    exit 1
fi


# Extract sample names
echo "Extracting sample names..."
bcftools query -l "$FILE_A" > samples_A.txt
bcftools query -l "$FILE_B" > samples_B.txt

# Generate sample name mapping
echo "Generating sample name mapping..."
awk '
    BEGIN {
        while (getline < "samples_A.txt") {
            original = $1;
            gsub(/^GM/, "NA", $1);
            if ($1 != original) {
                print original "\t" $1;
            }
            gsub(/^NA/, "GM", original);
            if ($1 != original) {
                print original "\t" $1;
            }
        }
    }
' > sample_name_mapping.txt

# Rename samples in File A
echo "Renaming samples in $FILE_A..."
bcftools reheader -s sample_name_mapping.txt "$FILE_A" > NW_flipped_50000_renamed.bcf

# Extract renamed sample names and find common samples
echo "Identifying common samples..."
bcftools query -l NW_flipped_50000_renamed.bcf > samples_A_renamed.txt
awk 'NR==FNR{samples[$1]=1; next} ($1 in samples)' samples_A_renamed.txt samples_B.txt > matching_samples.txt

# Isolate relevant samples
echo "Isolating relevant samples..."
bcftools view -S matching_samples.txt NW_flipped_50000_renamed.bcf -O b > A_restricted.bcf
bcftools view -S matching_samples.txt "$FILE_B" -O b > B_restricted.bcf

# Index the BCF files
echo "Indexing..."
tabix -p bcf A_restricted.bcf
sleep 5
tabix -p bcf B_restricted.bcf

# Concatenate the BCF files
echo "Concatenating files..."
bcftools concat A_restricted.bcf B_restricted.bcf -O z -a -o "$OUTPUT"

# Cleanup intermediate files
echo "Cleaning up..."
rm samples_A.txt samples_B.txt sample_name_mapping.txt NW_flipped_50000_renamed.bcf samples_A_renamed.txt matching_samples.txt A_restricted.bcf A_restricted.bcf.csi B_restricted.bcf B_restricted.bcf.csi

echo "All done! Merged file is saved as $OUTPUT."

