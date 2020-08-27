#!/usr/bin/awk -f
#usage: zcat file.fastq.gz | awk -v demuxmap=index_map.tsv -v odir=output_directory -f fastq_demux.awk -

BEGIN{
  FS=OFS="\t"; # set output separator
  while(getline < demuxmap) map[$1]=$2  # Load demultiplex sequences mapping
}

# Store 4 lines of a read into array r[]
{r[NR%4] = $0}

# Every time 4 lines are read
(NR%4==0) {
  # check FASTQ format
  if (r[1] !~ /^@/) exit 1;
  if (r[3] !~ /^[+]/) exit 2;

  #umi = substr(r[2],95,12);
  #barcode = substr(r[2],1,32);
  idx = substr(r[2],89,6);
  if (idx in map) {of = map[idx] ".fastq.gz"} else {of = "UNMAP.fastq.gz"}
  
  print r[1] "\n" r[2] "\n" r[3] "\n" r[0] | "pigz > " odir "/" of
  n[of]++
}

END {
  for(i in n) print i,n[i] > odir "/demux_report.txt"
}
