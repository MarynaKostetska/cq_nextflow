//General workflow

nextflow.enable.dsl=2

params.out = "${launchDir}/output"

// Download the Fasts file with a link
process downloadFile {
    publishDir params.out, mode: "copy", overwrite: true
    output:
        path "batch1.fasta"
    """
    wget https://tinyurl.com/cqbatch1 -O batch1.fasta
    """
}

// Count a general number of sequences and put the result in numseqs.txt file
process countSequences {
    publishDir params.out, mode: 'copy', overwrite: true
    input:
        path infile
    output:
        path "numseq*"
"""
grep "^>" ${infile} | wc -l > numseqs.txt
"""
}

// Split siquences from the fasta to differrent files with name split_0
process splitSequences {
  publishDir params.out, mode: "copy", overwrite: true
  input:
    path infile 
  output:
    path "seq_*.fasta"
  """
  split -d -l 2 $infile seq_
  for file in seq_*; do mv "\$file" "\$file.fasta"; done
  """
}
// Count a number of bases in each sequence
process countBases {
    input:
        infile
    output:
        path "${infile.getSimpleName()}.basescount"
"""
grep -v "^>" ${infile} | wc -m ${infile.getSimpleName()}.basescount
"""
}

// Count a number of "GCCGCG" repeats in each sequence (треба додати коментарі по коротких командах)
process countRepeats {
	publishDir params.out, mode: "copy", overwrite: true
	input:
		path infile
	output:
		path "${infile.getSimpleName()}.repeatcounts"
"""
echo -n "${infile.getSimpleName()}" | cut -d "_" -f 2 > ${infile.getSimpleName()}.repeatcounts
echo -n ", " >> ${infile.getSimpleName()}.repeatcounts
grep -o "GCCGCG" ${infile} | wc -l >> ${infile.getSimpleName()}.repeatcounts
"""
}

// Create a report
process makeReport {
    publishDir params.out, mode: "copy", overwrite: true
    input:
        path infile
    output:
        path "finalcount.csv"
"""
cat * > count.csv
echo "# Sequence number, repeats" > finalcount.csv
cat count.csv >> finalcount.csv
"""
}

workflow {
  downloadFile | splitSequences | flatten | countRepeats | collect | makeReport
}

/*
$ nextflow run <ім'я_вашого_скрипту>.nf
$ cd output
$ cat batch1.fasta
$ ls split_*
$ cat split_0
$ cd output
*/

