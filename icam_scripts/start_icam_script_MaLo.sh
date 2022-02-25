

conf=/flash/home/chritzel/nextflow/vcf2neofox/icam_scripts/ntpconf

 
function run {
python /flash/home/chritzel/nextflow/vcf2neofox/icam_scripts/nucleotide_to_peptide.py \
                                                              --mhc_I_selection=$m1 \
                                                              --mhc_II_selection=$m2 \
                                                              --gene_expression=$g \
                                                              --exon_expression="" \
                                                              --mutation_set_output=$o \
                                                              $conf \
                                                              $vcf vcf1
python /flash/home/chritzel/nextflow/vcf2neofox/icam_scripts/peptide_filter.py $o.transcript $o.transcript.squished "" 500000
python /flash/home/chritzel/nextflow/vcf2neofox/icam_scripts/calc_allele_distribution.py \
    --jar  /flash/home/chritzel/nextflow/vcf2neofox/icam_scripts/ReadCounts.jar \
    --java /code/iCaM2/SOUP/jdk1.7.0/bin/java --tmpfile tmp_"$id"_tmp --infile_has_header --add_to_file --prefix RNA_ --samples r $o $rnabam > $o.freq

}

 

 

patientid=test123
vcf=/scratch/info/projects/MASTER_Mainz/reports/MT000014_switched/MIRACUM/WES/somaticGermline_MT000014_VC.output.snp.Somatic.hc.vcf
g=/scratch/info/projects/CM01_iVAC/RB_0004/RB0401_01_001/scratch/scratch/iCaM_output/out.genes.norm.tab
rnabam=/scratch/info/projects/MASTER_Mainz/reports/MT000014_switched/MIRACUM/WES/somaticGermline_MT000014_TD_output.sort.filtered.rmdup.realigned.fixed.recal.bam
o="$patientid"_test

#HLA typen
m1="" # A,B
m2="" # DRB,DQBq

run

