#!/bin/bash


CHROM_SIZES="mm10.chrom.sizes"

declare -A groups=(
  ["0h"]="A0.bw B0.bw C0.bw"
  ["1K"]="A1K.bw B1K.bw C1K.bw"
  ["1M"]="A1M.bw B1M.bw C1M.bw"
  ["6K"]="A6K.bw B6K.bw C6K.bw"
  ["6M"]="A6M.bw B6M.bw C6M.bw"
  ["24K"]="A24K.bw B24K.bw C24K.bw"
  ["24M"]="A24M.bw B24M.bw C24M.bw"
)

for group in "${!groups[@]}"; do
  echo "Averaging ${group}..."
  wiggletools mean ${groups[$group]} > ${group}_mean.wig
  wigToBigWig ${group}_mean.wig $CHROM_SIZES ${group}_mean.bw
  rm ${group}_mean.wig  # Remove the intermediate WIG file
  echo "Done: ${group}_mean.bw"
done

echo "Averaged bigWigs written."
