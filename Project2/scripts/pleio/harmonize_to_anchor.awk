BEGIN { FS=OFS="\t" }
NR==FNR { a1[$1]=toupper($2); a2[$1]=toupper($3); next }
FNR==1 { next }
{
  snp=$1; A1=toupper($2); A2=toupper($3); Z=$4; N=$5
  if (!(snp in a1)) next

  r1=a1[snp]; r2=a2[snp]
  amb=((A1=="A"&&A2=="T")||(A1=="T"&&A2=="A")||(A1=="C"&&A2=="G")||(A1=="G"&&A2=="C"))
  if (amb) next

  if (A1==r1 && A2==r2) {
    print snp,r1,r2,Z,N
  } else if (A1==r2 && A2==r1) {
    print snp,r1,r2,-Z,N
  }
}
