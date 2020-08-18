!#/usr/bin/env nextflow

// enable dsl 2
nextflow.enable.dsl = 2
/*
* Pipeline Input Parameter
*/
params.reads = "test-datasets/testdata/*{1,2}.fastq.gz"
Channel.fromFilePairs(params.reads).into{reads_ch; reads_fastqc_ch}
params.fasta = "test-datasets/reference/genome.fa"
params.gtf = "test-datasets/reference/genome.gtf"
params.outdir = "results"
params.aligner = "bwa"
params.minNonOverlap = 10




name = "$workflow.runName"



log.info """\
      e c c D N A - N F  P I P E L I N E
      ==================================
      Name:   ${name}
      Genome: ${params.fasta}
      Reads:  ${params.reads}
      Outdir: ${params.outdir}
      ==================================
      """
      .stripIndent()

Channel.fromPath(params.fasta).into{fasta_for_bwaindex_ch; fasta_ch}
Channel.fromPath(params.gtf).set{bwa_index_gtf}

fasta_ch.view()

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
        file zipped_fasta

        output:
        file "*.{fa,fn,fna,fasta}" into ch_fasta_for_bwaindex

        script:
        rm_zip = zipped_fasta - '.gz'
        """
        pigz -f -d -p ${task.cpus} $zipped_fasta
        """
        }
    } else {
    fasta_for_indexing = Channel
    .fromPath("${params.fasta}", checkIfExists: true)
    .set{ ch_fasta_for_bwaindex }

    lastPath = params.fasta.lastIndexOf(File.separator)
    bwa_base = params.fasta.substring(lastPath+1)
}


process fastqc {
  tag "${sample_id}"

  publishDir "results/fastqc", mode: "copy"

  input:
  tuple val(sample_id), file(reads_file) from reads_fastqc_ch

  output:
  file "fastqc_${sample_id}_logs" into fastqc_ch

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q $reads_file
  """
}

process makeBWAindex {
  // cpus 2 <- use 2 cpus 
  tag "${fasta}"
  //publishDir "results/bwaIndex", mode: "copy"
  input:
  path fasta

  output:
  file "BWAindex" 

  script:
  """
  bwa index -a bwtsw $fasta
  mkdir BWAindex && mv ${fasta}* BWAindex
  """
}

// $task.cpus -> how many cpus used

process bwamem {
  publishDir  "results/mapping/bwamem", mode: "copy"
  tag "${sample_id}"
  when: params.aligner == "bwa"

  input:
  tuple val(sample_id), file(reads_file)
  // file index from bwa_index_ch.collect()
  file index

  output:
  file "*" into ch_output_from_bwamem
  val (sample_id) into sample_id_ch2
  //file "*.{bai,csi}" into ch_outputindex_from_bwamem

  script:
  index_path = "${index}/${bwa_base}"
  """
    bwa mem $index_path $reads_file > "${sample_id}".mapped.sam
  """
}

process samblaster {
  publishDir  "results/samblaster", mode: "copy"
  tag "${sample_id}"

  input:
  file mapped_reads from ch_output_from_bwamem
  val sample_id from sample_id_ch2

  output:
  tuple val(sample_id), file("${sample_id}.disc.bam"), file("${sample_id}.split.bam"), file("${sample_id}.concordant.bam") into samblaster_bam_ch
//  file "*.disc.sam" into samblaster_disc_ch
//  file "*.split.sam" into samblaster_split_ch
//  file "*.unmap.fastq" into samblaster_unmap_ch
//  file "${sample_id}.sam" into samblaster_sam_ch
//  val (sample_id) into sample_id_ch3

  script:
  """
  samtools view -h $mapped_reads | \
  samblaster -e --minNonOverlap $params.minNonOverlap -d "${sample_id}.disc.sam" -s "${sample_id}.split.sam" -u "${sample_id}.unmap.fastq" > "${sample_id}.sam"

  samtools view -bS "${sample_id}.sam" -o "${sample_id}.bam"
  samtools sort -O bam -o "${sample_id}.sorted.bam" "${sample_id}.bam"
  samtools index "${sample_id}.sorted.bam"

  samtools view -bS ${sample_id}.disc.sam > "${sample_id}.disc.bam"
  samtools view -bS ${sample_id}.split.sam > "${sample_id}.split.bam"
  samtools view -hf 0x2 ${sample_id}.sorted.bam -bS > "${sample_id}.concordant.bam"
  """
}

//process samtoolsIndex {
//  publishDir "results/samtools/index", mode: "copy"
//  tag "${sample_id}"
//
//  input:
//  file sam from samblaster_sam_ch
//  val sample_id from sample_id_ch3
//
//  output:
//  file "${sample_id}.bam" into bam_file_ch
//  file "*.sorted.bam" into bam_sorted_ch
//  file "*.sorted.bam.bai" into bam_sorted_index_ch
//  val sample_id into sample_id_ch4
//
//  script:
//  """
//  samtools view -bS $sam -o "${sample_id}.bam"
//  samtools sort -O bam -o "${sample_id}.sorted.bam" "${sample_id}.bam"
//  samtools index "${sample_id}.sorted.bam"
//  """
//}

//process samToBam {
//  publishDir "results/samtools/bamFiles", mode: "copy"
//  tag "${sample_id}"
//
//  input:
//  file disc from samblaster_disc_ch
//  file split from samblaster_split_ch
//  file bam_sorted from bam_sorted_ch
//  val sample_id from sample_id_ch4
//
//  output:
//  file "*.disc.bam" into disc_bam_ch
//  file "*.split.bam" into split_bam_ch
//  file "*concordant.bam" into concordant_bam_ch
//  val sample_id into sample_id_ch5
//  //file "*.unmap.bam" into unmap_bam_ch
//
//  script:
//  """
//  samtools view -bS $disc > "${sample_id}.disc.bam"
//  samtools view -bS $split > "${sample_id}.split.bam"
//  samtools view -hf 0x2 $bam_sorted -bS > "${sample_id}.concordant.bam"
//  """
//}

//process bamToBed {
//  publishDir "results/bedtools/bedFiles", mode: "copy"
//  tag "${sample_id}"
//  echo true
//
//  input:
//
//  output:
//
//  script:
//  """
//  """
//}

process bam_to_bed {
  publishDir "results/bedFiles", mode: "copy"

  input:
  //file split_txt from split_txt_ch
  //file concordant_txt from concordant_txt_ch
  tuple val(sample_id), file(disc_bam), file(split_bam), file(concordant_bam) from samblaster_bam_ch

  output:
  tuple val(sample_id),
        file("${sample_id}.split.txt"),
        file("${sample_id}.concordant.txt") into split_concordant_ch

  file "*.disc.txt" into disc_txt_ch

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


//  awk '{print \$4}' ${sample_id}.split.txt | sort | uniq -c > "${sample_id}.split.id-freq.txt"
//  awk '\$1=="2" {print \$2}' "${sample_id}.split.id-freq.txt" > "${sample_id}.split.id-freq2.txt"
//  awk '\$1=="4" {print \$2}' "${sample_id}.split.id-freq.txt" > "${sample_id}.split.id-freq4.txt"
//
//  awk '{print \$4}' ${sample_id}.concordant.txt | sort | uniq -c > "${sample_id}.concordant.id-freq.txt"
//  awk '\$1=="3" {print \$2}' "${sample_id}.concordant.id-freq.txt" > "${sample_id}.concordant.id-freq3.txt"
//  awk '\$1>3 {print \$2}' "${sample_id}.concordant.id-freq.txt" > "${sample_id}.concordant.id-freqGr3.txt"
process circle_finder {
  publishDir "results/circle_finder", mode: "copy"

  input:
  tuple val(sample_id),
        file(split_txt),
        file(concordant_txt) from split_concordant_ch

  output:
  file "*" into circle_finder_ch

  //when:
  //split_id_freq2.size() > 0 & split_id_freq4.size() > 0 & concordant_id_freq3.size() > 0 & concordant_id_freqGr3.size() > 0

  script:
  """
  circle_finder.sh $sample_id \
                   $split_txt \
                   $concordant_txt
  """
}

//  file "**.split_freq2.txt" into split_freq2_ch
//  file "**.split_freq4.txt" into split_freq4_ch
//
//  file "**.concordant_freq3.txt" into concordant_freq3_ch
//  file "**.concordant_freqGr3.txt" into concordant_freqGr3_ch
//  file "**.split_freq2.oneline.txt" into split_freq2_oneline_ch
//  file "**.split_freq4.oneline.txt" into split_freq4_oneline_ch
//  file "**.concordant_freq3.2SPLIT-1M.txt" into conc_freq3_2SPLIT_1M_ch
//  file "**.concordant_freq3.2SPLIT-1M.inoneline.txt" into conc_freq3_2SPLIT_1M_oneline_ch
//  file "**.microDNA-JT.txt" into microDNA_JT_ch
process multiqc {
  publishDir "results/multiqc", mode: "copy"

  input:
  path "*" from fastqc_ch.collect()

  output:
  path 'multiqc_report.html'

  script:
  """
  multiqc .
  """
}

workflow.onComplete {
  log.info ( workflow.success ? "\n Done! Open the multiqc report in your browser --> $params.outdir/multiqc/multiqc_report.html\n" : 
                                "\n Error. Pipeline execution stopped with the following message: ${workflow.errorMessage}\n")
}

// define dsl2 variables
workflow {
  fastqc(read_pairs_ch)
  makeBWAindex(params.fasta)
  bwamem(read_pairs_ch, makeBWAindex.out.collect())
}

//  mkdir "${sample_id}"
//  grep -w -Ff $split_id_freq2 $split_txt > "${sample_id}/${sample_id}.split_freq2.txt"
//  grep -w -Ff $split_id_freq4 $split_txt > "${sample_id}/${sample_id}.split_freq4.txt"
//
//  grep -w -Ff $concordant_id_freq3 $concordant_txt > "${sample_id}/${sample_id}.concordant_freq3.txt"
//  grep -w -Ff $concordant_id_freqGr3 $concordant_txt > "${sample_id}/${sample_id}.concordant_freqGr3.txt"
//
//  sed 'N;s/\\n/\t/' "${sample_id}/${sample_id}.split_freq2.txt" > "${sample_id}/${sample_id}.split_freq2.oneline.txt"
//  sed 'N;s/\\n/\t/' "${sample_id}/${sample_id}.split_freq4.txt" > "${sample_id}/${sample_id}.split_freq4.oneline.txt"
//
//  awk '\$1==\$10 && \$7==\$16 && \$6>0 && \$15>0 {print \$4} ' "${sample_id}/${sample_id}.split_freq2.oneline.txt" > \
//    "${sample_id}/${sample_id}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt"
//
//  grep -w -Ff "${sample_id}/${sample_id}.split_freq2.oneline.S-R-S-CHR-S-ST.ID.txt" "${sample_id}/${sample_id}.concordant_freq3.txt" > \
//    "${sample_id}/${sample_id}.concordant_freq3.2SPLIT-1M.txt"
//
//  awk 'BEGIN{FS=OFS="\t"} {gsub("M", " M ", \$8)} 1' "${sample_id}/${sample_id}.concordant_freq3.2SPLIT-1M.txt" | \
//    awk 'BEGIN{FS=OFS="\t"} {gsub("S", " S ", \$8)} 1' | \
//    awk 'BEGIN{FS=OFS="\t"} {gsub("H", " H ", \$8)} 1' | \
//    awk 'BEGIN{FS=OFS="\t"} {gsub("S", " S ", \$8)} 1' | \
//    awk 'BEGIN{FS=OFS="\t"} {gsub("H", " H ", \$8)} 1' | \
//    awk 'BEGIN{FS=OFS=" "} {if ((\$9=="M" && \$NF=="H") || \
//      (\$9=="M" && \$NF=="S"))  {printf ("%s\tfirst\\n",\$0)} else if ((\$9=="S" && \$NF=="M") || \
//      (\$9=="H" && \$NF=="M")) {printf ("%s\tsecond\\n",\$0)} else  {printf ("%s\tconfusing\\n",\$0)}}' | \
//    awk 'BEGIN{FS=OFS="\t"} {gsub(" ", "", \$8)} 1' | \
//    awk '{printf ("%s\t%d\\n",\$0,(\$3-\$2)+1)}' | sort -k4,4 -k10,10n | \
//    sed 'N;N;s/\\n/\\t/g' | awk '{if (\$5 == \$15) {print \$0} else if (( \$5 == "1" && \$15 == "2" && \$25 == 1) || \
//      (\$5 == "2" && \$15 == "1" && \$25 == "2")) \
//      {printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\\n", \
//      \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$21,\$22,\$23,\$24,\$25,\$26,\$27,\$28,\$29,\$30,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20)} \
//      else if ((\$5=="1" && \$15=="2" && \$25=="2") || (\$5=="2" && \$15=="1" && \$25=="1")) \
//      {printf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%d\\n", \
//      \$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$21,\$22,\$23,\$24,\$25,\$26,\$27,\$28,\$29,\$30,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10)}}' \
//      > "${sample_id}/${sample_id}.concordant_freq3.2SPLIT-1M.inoneline.txt"
//
//  awk '\$1==\$11 && \$1==\$21 && \$7==\$17'  "${sample_id}/${sample_id}.concordant_freq3.2SPLIT-1M.inoneline.txt" | \
//    awk '(\$7=="+" && \$27=="-") || (\$7=="-" && \$27=="+")' | \
//    awk '{if (\$17=="+" && \$19=="second" && \$12<\$2 && \$22>=\$12 && \$23<=\$3) {printf ("%s\t%d\t%d\n",\$1,\$12,\$3)} else if \
//      (\$7=="+" && \$9=="second" && \$2<\$12 && \$22>=\$2 && \$23<=\$13) {printf ("%s\t%d\t%d\\n",\$1,\$2,\$13)} else if \
//      (\$17=="-" && \$19=="second" && \$12<\$2 && \$22>=\$12 && \$23<=\$3) {printf ("%s\t%d\t%d\\n",\$1,\$12,\$3)} else if \
//      (\$7=="-" && \$9=="second" && \$2<\$12 && \$22>=\$2 && \$23<=\$13) {printf ("%s\t%d\t%d\\n",\$1,\$2,\$13)} }' | \
//    sort | uniq -c | awk '{printf ("%s\t%d\t%d\t%d\\n",\$2,\$3,\$4,\$1)}' > "${sample_id}/${sample_id}.microDNA-JT.txt"