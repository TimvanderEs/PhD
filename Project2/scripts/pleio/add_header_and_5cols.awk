BEGIN { FS=OFS="\t" }
NR==1 {
  if ($1=="SNP" && $2=="A1" && $3=="A2") next
}
{ print $1,$2,$3,$4,$5 }
