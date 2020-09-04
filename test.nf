#!/usr/bin/env nextflow

// enable dsl 2
nextflow.enable.dsl = 2



log.info """\
      e c c D N A - N F  P I P E L I N E
      ==================================
      Name:   ${workflow.runName}
      Genome: ${params.fasta}
      Reads:  ${params.reads}
      Outdir: ${params.outdir}
      Save Trimmed: ${params.saveTrimmed}
      ==================================
      """
      .stripIndent()


if (!params.bwa_index) {
Channel.fromPath(params.fasta).set{fasta_ch}
if ( params.fasta.isEmpty () ){
    exit 1, "Please specify --fasta with the path to your reference"
} else if("${params.fasta}".endsWith(".gz")){
    //Put the zip into a channel, then unzip it and forward to downstream processes. DONT unzip in all steps, this is inefficient as NXF links the files anyways from work to work dir
    zipped_fasta = file("${params.fasta}")

    rm_gz = params.fasta - '.gz'
    lastPath = rm_gz.lastIndexOf(File.separator)
    bwa_base = rm_gz.substring(lastPath+1)

    process unzip_reference{
        tag "${zipped_fasta}"

        input:
        file(zipped_fasta)

        output:
        file("*.{fa,fn,fna,fasta}") into fasta_for_bwaindex_ch

        script:
        rm_zip = zipped_fasta - '.gz'
        """
        pigz -f -d -p ${task.cpus} $zipped_fasta
        """
        }
    } else {
    fasta_for_indexing = Channel
    .fromPath("${params.fasta}", checkIfExists: true)
    .set{ fasta_for_bwaindex_ch }

    lastPath = params.fasta.lastIndexOf(File.separator)
    bwa_base = params.fasta.substring(lastPath+1)
  } 
} else {
    lastPath = params.bwa_index.lastIndexOf(File.separator) 
    bwa_dir  = params.bwa_index.substring(0,lastPath+1)
    bwa_base = params.bwa_index.substring(lastPath+1)
    Channel .fromPath(bwa_dir, checkIfExists: true) .set { ch_bwa_index }
}


if (!params.skip_fastqc){
  process fastqc {
    tag "${sample_id}"
    label "low_memory"

    publishDir "${params.outdir}/fastqc", mode: "copy"

    input:
    tuple val(sample_id), file(reads)

    output:
    file("${sample_id}_logs")

    script:
    """
    mkdir ${sample_id}_logs
    fastqc --quiet --threads $task.cpus -o ${sample_id}_logs -q $reads
    """
  }
}

process trim_galore {
  label 'low_memory'
  tag "$sample_id"
  publishDir "${params.outdir}/trim_galore", mode: 'copy',
   saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else if (params.saveTrimmed) filename
                else if (!params.saveTrimmed) null
   }

  input:
  tuple val(sample_id), file(reads)

  output:
  tuple val(sample_id), file("*fq.gz") 
  file "*trimming_report.txt"
  file "*_fastqc.{zip,html}"

  script:
  c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
  c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
  tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
  tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
  nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''

  """
  trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $nextseq $reads
  """
}

process makeBWAindex {
  tag "${fasta}"
  label "high_memory"

  input:
  path(fasta)

  output:
  file("BWAindex")

  script:
  """
  bwa index -a bwtsw $fasta
  mkdir BWAindex && mv ${fasta}* BWAindex
  """
}

// $task.cpus -> how many cpus used

process bwamem {
  //publishDir "${params.outdir}/bwamem/${sample_id}"
  tag "${sample_id}"
  label "high_memory"
  when: params.aligner == "bwa"

  input:
  tuple val(sample_id), file(reads)
  file(index)

  output:
  tuple val(sample_id), file("*.mapped.sam")

  script:
  index_path = "${index}/${bwa_base}"
  """
    bwa mem -t $task.cpus $index_path $reads > "${sample_id}".mapped.sam
  """
}

process samblaster {
  publishDir  "${params.outdir}/samblaster/${sample_id}", mode: "copy"
  label "mid_memory"
  tag "${sample_id}"

  input:
  tuple val(sample_id), file(mapped_reads)

  output:
    tuple val(sample_id),
          file("${sample_id}.disc.bam"),
          file("${sample_id}.split.bam"),
          file("${sample_id}.concordant.bam")
    file("${sample_id}.sorted.bam")
    file("${sample_id}.sorted.bam.bai")

  script:
  """
  samtools view --threads $task.cpus -h $mapped_reads | \
  samblaster -e --minNonOverlap $params.minNonOverlap -d "${sample_id}.disc.sam" -s "${sample_id}.split.sam" -u "${sample_id}.unmap.fastq" > "${sample_id}.sam"

  samtools view --threads $task.cpus \
    -bS "${sample_id}.sam" -o "${sample_id}.bam"
  samtools sort --threads $task.cpus -O bam -o "${sample_id}.sorted.bam" "${sample_id}.bam"
  samtools index -@ $task.cpus "${sample_id}.sorted.bam"

  samtools view --threads $task.cpus -bS ${sample_id}.disc.sam > "${sample_id}.disc.bam"
  samtools view --threads $task.cpus -bS ${sample_id}.split.sam > "${sample_id}.split.bam"
  samtools view --threads $task.cpus -hf 0x2 ${sample_id}.sorted.bam -bS > "${sample_id}.concordant.bam"
  """
}

process bam_to_bed {
  publishDir "${params.outdir}/bedFiles/${sample_id}", mode: "copy"
  tag "${sample_id}"
  label "mid_memory"

  input:
  tuple val(sample_id),
        file(disc_bam),
        file(split_bam),
        file(concordant_bam)

  output:
  tuple val(sample_id),
        file("${sample_id}.split.txt"),
        file("${sample_id}.concordant.txt")

  script:
  """
  bedtools bamtobed -cigar -i ${split_bam} | sed -e 's/_2\\/2/ 2/g' | \
    sed -e 's/_1\\/1/ 1/g' |
    awk '{printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\\n", \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8}' |
    awk 'BEGIN{FS=OFS="\t"} {gsub("M", " M ", \$8)} 1' | \
    awk 'BEGIN{FS=OFS="\t"} {gsub("S", " S ", \$8)} 1' | \
    awk 'BEGIN{FS=OFS="\t"} {gsub("H", " H ", \$8)} 1' | \
    awk 'BEGIN{FS=OFS=" "} {if ((\$9=="M" && \$NF=="H") || \
      (\$9=="M" && \$NF=="S"))  {printf ("%s\tfirst\\n",\$0)} else if ((\$9=="S" && \$NF=="M") || \
      (\$9=="H" && \$NF=="M")) {printf ("%s\tsecond\\n",\$0)} }' | \
    awk 'BEGIN{FS=OFS="\t"} {gsub(" ", "", \$8)} 1' > '${sample_id}.split.txt'

  bedtools bamtobed -cigar -i ${concordant_bam} | \
    sed -e 's/\\// /g' | \
    awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8)}' > '${sample_id}.concordant.txt'

  bedtools bamtobed -cigar -i ${disc_bam} | \
    sed -e 's/\\// /g' | \
    awk '{printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8)}' > '${sample_id}.disc.txt'
  """
}

process circle_finder {
  publishDir "${params.outdir}/circle_finder/${sample_id}", mode: "copy"

  echo true
  tag "${sample_id}"
  label "mid_memory"

  input:
  tuple val(sample_id),
        file(split_txt),
        file(concordant_txt)

  shell:
  '''
  #!/usr/bin/env bash
  set -u

  sample_id=!{sample_id}
  split=!{split_txt}
  concordant=!{concordant_txt}

  # Function to output an error if files do not exist
  file_exists () {
      [[ ! -s $1 ]] && \
      echo "
  ERROR - CIRCLE_FINDER - $sample_id
  ===================================
  Error $1 does not exist or is empty.
  Stopped circle_finder.
  No circular DNA was identified.
  ===================================
  " && exit
  }

  awk '{print $4}' ${split} | sort | uniq -c > ${sample_id}.split.id-freq.txt
  #This file "${sample_id}.split.id-freq.txt" will be used for collecting split id that have frequency equal to 4.
  awk '$1=="2" {print $2}' ${sample_id}.split.id-freq.txt > ${sample_id}.split.id-freq2.txt
  awk '$1=="4" {print $2}' ${sample_id}.split.id-freq.txt > ${sample_id}.split.id-freq4.txt

  awk '{print $4}' ${concordant} | sort | uniq -c > ${sample_id}.concordant.id-freq.txt
  #The following command will chose (may not be always true) one concordant and 2 split read

  awk '$1=="3" {print $2}' ${sample_id}.concordant.id-freq.txt > ${sample_id}.concordant.id-freq3.txt
  awk '$1>3 {print $2}' ${sample_id}.concordant.id-freq.txt > ${sample_id}.concordant.id-freqGr3.txt

  # check if output files exists and are not empty
  file_exists ${sample_id}.concordant.id-freq3.txt
  file_exists ${sample_id}.concordant.id-freqGr3.txt

  grep -w -Ff ${sample_id}.split.id-freq2.txt ${split} > ${sample_id}.split_freq2.txt
  grep -w -Ff ${sample_id}.split.id-freq4.txt ${split} > ${sample_id}.split_freq4.txt

  # check if output files exist and are not empty
  file_exists ${sample_id}.split_freq2.txt
  file_exists ${sample_id}.split_freq4.txt

  #Selecting concordant pairs that were 1) mapped uniquely and 2) mapped on more than one loci (file "freqGr3.txt")
  grep -w -Ff ${sample_id}.concordant.id-freq3.txt ${concordant} > ${sample_id}.concordant_freq3.txt
  grep -w -Ff ${sample_id}.concordant.id-freqGr3.txt ${concordant} > ${sample_id}.concordant_freqGr3.txt

  # check if output files exist and are not empty
  file_exists ${sample_id}.concordant_freq3.txt
  file_exists ${sample_id}.concordant_freqGr3.txt

  #Step 7: Putting split read with same id in one line
  sed 'N;s/\\n/\\t/' ${sample_id}.split_freq2.txt > ${sample_id}.split_freq2.oneline.txt
  sed 'N;s/\\n/\\t/' ${sample_id}.split_freq4.txt > ${sample_id}.split_freq4.oneline.txt

  # check if output files exist and are not empty
  file_exists ${sample_id}.split_freq2.oneline.txt
  file_exists ${sample_id}.split_freq4.oneline.txt

  #Step 8: Split reads map on same chromosome and map on same strand. Finally extracting id (split read same chromosome, split read same strand), collecting all the split reads that had quality >0
  awk '$1==$10 && $7==$16 && $6>0 && $15>0 {print $4} ' ${sample_id}.split_freq2.oneline.txt > \
      ${sample_id}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt

  # check if output files exist and are not empty
  file_exists ${sample_id}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt

  #Step 9: Based on unique id I am extracting one continuously mapped reads and their partner mapped as split read (3 lines for each id) 
  grep -w -Ff "${sample_id}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt" "${sample_id}.concordant_freq3.txt" > \
      "${sample_id}.concordant_freq3.2SPLIT-1M.txt"

  # check if output files exist and are not empty
  file_exists ${sample_id}.concordant_freq3.2SPLIT-1M.txt

  #Step 10: Sorting based on read-id followed by length of mapped reads.
  awk 'BEGIN{FS=OFS="\\t"} {gsub("M", " M ", $8)} 1' ${sample_id}.concordant_freq3.2SPLIT-1M.txt | \
      awk 'BEGIN{FS=OFS="\\t"} {gsub("S", " S ", $8)} 1' | \
      awk 'BEGIN{FS=OFS="\\t"} {gsub("H", " H ", $8)} 1' | \
      awk 'BEGIN{FS=OFS=" "} {if (($9=="M" && $NF=="H") || \
          ($9=="M" && $NF=="S"))  {printf ("%s\\tfirst\\n",$0)} \
          else if (($9=="S" && $NF=="M") || ($9=="H" && $NF=="M")) {printf ("%s\\tsecond\\n",$0)} \
          else  {printf ("%s\\tconfusing\\n",$0)}}' | \
      awk 'BEGIN{FS=OFS="\\t"} {gsub(" ", "", $8)} 1' | \
      awk '{printf ("%s\\t%d\\n",$0,($3-$2)+1)}' | \
      sort -k4,4 -k10,10n | sed 'N;N;s/\\n/\\t/g' | \
      awk '{if ($5==$15) {print $0}  \
          else if (($5=="1" && $15=="2" && $25=="1") || ($5=="2" && $15=="1" && $25=="2")) \
              {printf ("%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20)} \
          else if (($5=="1" && $15=="2" && $25=="2") || ($5=="2" && $15=="1" && $25=="1")) \
          {printf ("%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\t%s\\t%d\\t%d\\t%s\\t%d\\t%d\\t%s\\t%s\\t%s\\t%d\\n", $11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10)} }' \
      > ${sample_id}.concordant_freq3.2SPLIT-1M.inoneline.txt

  # check if output files exist and are not empty
  file_exists ${sample_id}.concordant_freq3.2SPLIT-1M.inoneline.txt

  #Step 11: Unique number of microDNA with number of split reads
  awk '$1==$11 && $1==$21 && $7==$17'  ${sample_id}.concordant_freq3.2SPLIT-1M.inoneline.txt | \
      awk '($7=="+" && $27=="-") || ($7=="-" && $27=="+")' | \
      awk '{if ($17=="+" && $19=="second" && $12<$2 && $22>=$12 && $23<=$3) {printf ("%s\\t%d\\t%d\\n",$1,$12,$3)} \
          else if ($7=="+" && $9=="second" && $2<$12 && $22>=$2 && $23<=$13) {printf ("%s\\t%d\\t%d\\n",$1,$2,$13)} \
          else if ($17=="-" && $19=="second" && $12<$2 && $22>=$12 && $23<=$3) {printf ("%s\\t%d\\t%d\\n",$1,$12,$3)} \
          else if ($7=="-" && $9=="second" && $2<$12 && $22>=$2 && $23<=$13) {printf ("%s\\t%d\\t%d\\n",$1,$2,$13)} }' | \
      sort | uniq -c | awk '{printf ("%s\\t%d\\t%d\\t%d\\n",$2,$3,$4,$1)}' > ${sample_id}.microDNA-JT.txt

  file_exists ${sample_id}.microDNA-JT.txt
  '''
}

process multiqc {
  publishDir "${params.outdir}/multiqc", mode: "copy"
  label "low_memory"

  input:
  path "*"

  output:
  path "${workflow.runName}_multiqc_report.html"

  script:
  """
  multiqc -n "${workflow.runName}_multiqc_report.html" .
  """
}

workflow.onComplete {
  log.info ( workflow.success ? "\n Done! Open the multiqc report in your browser --> $params.outdir/multiqc/${workflow.runName}_multiqc_report.html\n" : 
                                "\n Error. Pipeline execution stopped with the following message: ${workflow.errorMessage}\n")
}

/*
* Pipeline Input Parameter
*/
Channel.fromFilePairs(params.reads).set{read_pairs_ch}
// define dsl2 variables
workflow {
  if (!params.bwa_index) {
    makeBWAindex(fasta_for_bwaindex_ch)
    ch_bwa_index = makeBWAindex.out.collect()
  }

  if (!params.skipTrimming && !params.skip_fastqc){
    fastqc(read_pairs_ch)
    trim_galore(read_pairs_ch)
    bwamem(trim_galore.out[0], ch_bwa_index.collect())
    multiqc(fastqc.out.mix(trim_galore.out[1]).collect())
  } else if (!params.skip_fastqc && params.skipTrimming) {
    fastqc(read_pairs_ch)
    bwamem(read_pairs_ch, ch_bwa_index.collect())
    multiqc(fastqc.out.collect())
  } else if (params.skip_fastqc && params.skipTrimming) {
    bwamem(read_pairs_ch, ch_bwa_index.collect())
  }

  samblaster(bwamem.out)
  bam_to_bed(samblaster.out[0])
  circle_finder(bam_to_bed.out[0])
}
