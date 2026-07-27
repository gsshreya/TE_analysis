#!/bin/bash

# combine all alignment files of samples from /gpfs/gidata/genomeindia/eric/bams, /gpfs/gidata/genomeindia/eric/nibmg_samples/bams and /gpfs/data/user/shweta_lab/data/bams/GI/markduplicates locations and CREATE SOFTLINKS at /gpfs/data/user/ruchitha/GI_bams_links/6k_samples.
# bash /gpfs/data/user/ruchitha/scripts/CreateSoftLinks.sh
bash /gpfs/data/user/ruchitha/scripts/CoverageCalculation/CreateSoftLinksForCoverageCalculation.sh

#markduplicates_dir=/gpfs/data/user/shweta_lab/data/bams/GI/markduplicates/
markduplicates_dir="/gpfs/data/user/ruchitha/GI_bams_links/7480_samples"
misc_files="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/misc_files/"
cd ${markduplicates_dir}
for sample in $(ls *_md.bam.bai | grep -v -E '_1_|_2_|_3_|cbr'); do if [ -e ${sample%.bai} ]; then echo ${sample}; fi; done > ${misc_files}markdups_samples.txt

#grep -F -f ${misc_files}Joint_call_samples.txt ${misc_files}markdups_samples.txt > ${misc_files}markdups_samples_final_callset.txt
grep -F -f /gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/BatchesRun/57_combined_torerunMELT_for7480samples.txt ${misc_files}markdups_samples.txt > ${misc_files}markdups_samples_final_callset.txt

Work_dir="/gpfs/data/user/shweta_lab/data/TE/TE_discovery/joint_call_samples/insertions/"
samples_file=${misc_files}markdups_samples_final_callset.txt

#create softlinks of new samples in the working directory
bash /gpfs/data/user/shweta_lab/data/TE/TE_discovery/scripts/create_softlinks_frommultiplelocations.sh

while true; do
  running_jobs=$(squeue -u ruchitha | wc -l) #change the username accordingly
  if [ $running_jobs -lt 4 ]; then #change the number of maximum jobs to run (num jobs + 1)
    echo "Sleeping for 20s before next submission"
    sleep 20s
    counter=0
    samples=""
    while read -r sample_name; do
      prefix=$(echo "$sample_name" | cut -d '.' -f 1)
      if [ -e "${Work_dir}${prefix}.bam.disc" ] || \
         ( [ -e "${Work_dir}${prefix}.LINE1.aligned.final.sorted.bam" ] && \
         [ -e "${Work_dir}${prefix}.SVA.aligned.final.sorted.bam" ] && \
         [ -e "${Work_dir}${prefix}.HERVK.aligned.final.sorted.bam" ] && \
         [ -e "${Work_dir}${prefix}.ALU.aligned.final.sorted.bam" ] ); then
        continue
      else
        samples+="${sample_name%.*} "
        ((counter++))
        if [ $counter -ge 19 ]; then  # number of samples per job. Max to 20.
          # Create the second script using sed
          cd /gpfs/data/user/shweta_lab/data/TE/TE_discovery/scripts/
          sed "s|{{SAMPLES}}|${samples}|g" /gpfs/data/user/ruchitha/MELT/CallingMELT/actualrun/scripts/MELT_template.sh > script_part2.sh
          # Submit the second script if there are fewer than 6 running jobs
          sbatch script_part2.sh
          rm -f script_part2.sh
          # Reset the counter and samples list
          counter=0
          samples=""
          sleep 10s
          break  # Break out of the inner loop
        fi
      fi
    done < "${samples_file}"
  else
    echo "Maximum number of jobs reached. Waiting..."
    sleep 1m  # Sleep for 1 minute before checking again
  fi
done