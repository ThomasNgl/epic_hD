FA=/mnt/thomas/matriline/results/GRCm39/mm39.fa
samtools faidx "$FA"

seq=$(samtools faidx "$FA" chrM | grep -v '^>' | tr -d '\n' | tr 'acgtn' 'ACGTN')

echo "$seq" | awk '{
  s = $0              # keep original
  t = $0              # copy for C counting

  c = gsub(/C/, "", t)  # count Cs by modifying the COPY, not s

  cp = 0
  for (i = 1; i < length(s); i++) {
    if (substr(s, i, 2) == "CG") cp++
  }

  print "chrM_C_bases=" c
  print "chrM_CpG_dinucs=" cp
  print "chrM_CpG_cytosines_both_strands=" 2*cp
}'
