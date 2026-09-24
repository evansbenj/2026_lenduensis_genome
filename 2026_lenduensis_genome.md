# Lenduensis assembly

Flye
```
file                       format  type  num_seqs        sum_len  min_len   avg_len    max_len     Q1      Q2      Q3  sum_gap      N50  Q20(%)  Q30(%)  GC(%)
lendu_flye_assembly.fasta  FASTA   DNA     76,089  5,798,399,534       77  76,205.5  3,319,733  4,471  15,894  62,834        0  323,978       0       0  39.64

```

HiFiasm
```
file                            format  type  num_seqs        sum_len  min_len   avg_len    max_len      Q1      Q2      Q3  sum_gap     N50  Q20(%)  Q30(%)  GC(%)
lendu_hifiasm.bp.hap1.p_ctg.fa  FASTA   DNA    120,309  6,104,902,235    1,372  50,743.5  1,122,594  21,231  36,482  63,539        0  71,031       0       0  39.73
```

Hifiasm is longer but with more scaffolds and lower N50.


# kmerz in windows

Make genomic windows:
```
cut -f1,2 lendu_hifiasm.bp.hap1.p_ctg.fa.fai > genome.sizes
WINDOW=100000
bedtools makewindows -g genome.sizes -w $WINDOW > windows.bed
bedtools getfasta -fi lendu_hifiasm.bp.hap1.p_ctg.fa -bed windows.bed -fo windows.fa
```
split the windows.fa into smaller filez
```
N=$(grep -c '^>' windows.fa)
PER=$(( (N + 49) / 50 ))
awk -v per="$PER" '
/^>/ {
rec++
file=sprintf("split_%02d.fa", int((rec-1)/per)+1)
}
{ print >> file }
' windows.fa
```
# make smaller windowz
```
awk -v N=50 '
BEGIN {
    file=1;
    perfile=7240;
}
/^>/ {
    if (++n > perfile && file < N) {
        file++;
        n=1;
    }
}
{
    print > sprintf("split_%02d.fa", file);
}
' input.fa
```

# Do intersectionz
```
#!/usr/bin/env bash
set -euo pipefail
INPUT_FA="${1:?Usage: $0 input.fa}"
QUERY_DB="in_lendu_nonopore_meryldb.out_but_not_SRR38285470_trim.R12_meryldb.out.meryl"
K=29
THREADS=${SLURM_CPUS_PER_TASK:-1}

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

process_record() {
    local fa="$1"
    local id="$2"

    rm -rf "$tmpdir/window.meryl" "$tmpdir/intersect.meryl"

    /home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl count \
        k=${K} \
        threads=${THREADS} \
        output "$tmpdir/window.meryl" \
        "$fa" >/dev/null

    /home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl intersect \
	threads=${THREADS} \
        "$tmpdir/window.meryl" \
        "$QUERY_DB" \
        output "$tmpdir/intersect.meryl" >/dev/null

    count=$(/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl statistics "$tmpdir/intersect.meryl" |
            awk '/distinct/ {print $2; exit}')

    printf "%s\t%s\n" "$id" "${count:-0}"
}

outfile="$tmpdir/current.fa"
id=""

while IFS= read -r line
do
    if [[ $line == ">"* ]]; then

        if [[ -n "$id" ]]; then
            process_record "$outfile" "$id"
        fi

        id=${line#>}
        printf "%s\n" "$line" > "$outfile"

    else
        printf "%s\n" "$line" >> "$outfile"
    fi

done < "$INPUT_FA" > "$INPUT_FA"_counts.tsv

if [[ -n "$id" ]]; then
    process_record "$outfile" "$id" >> "$INPUT_FA"_counts.tsv
fi
```

# Get stats
```
awk '
{
    split($1,a,":")
    split(a[2],b,"-")
    len=b[2]-b[1]
    x=$NF/len

    v[n++] = x
    sum += x

    if (n==1 || x<min) min=x
    if (n==1 || x>max) max=x
}
END {
    mean = sum/n

    ss = 0
    for (i=0; i<n; i++)
        ss += (v[i]-mean)^2

    sd = sqrt(ss/(n-1))
    se = sd/sqrt(n)

    ci_low  = mean - 1.96*se
    ci_high = mean + 1.96*se

    printf "N\t%d\n", n
    printf "mean\t%.10f\n", mean
    printf "min\t%.10f\n", min
    printf "max\t%.10f\n", max
    printf "sd\t%.10f\n", sd
    printf "95%%CI\t%.10f\t%.10f\n", ci_low, ci_high
}' file.txt
```
# Concatenate countz
```
cat *counts.txt > all_countz.txt
```

# Summarize counts:
```
awk '{
  split($1,a,":");
  split(a[2],b,"-");
  len=b[2]-b[1];
  print $0, $2/len
}' all_countz.txt | sort -k3,3nr | head -5
```
