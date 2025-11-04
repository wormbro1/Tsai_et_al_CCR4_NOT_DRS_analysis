#!/bin/bash
#BSUB -J dorado_rep1
#BSUB -o /beevol/home/ouyangjo/logs/z_std_%J.out
#BSUB -e /beevol/home/ouyangjo/logs/z_std_%J.err
#BSUB -R "rusage[mem=64] span[hosts=1]"
#BSUB -gpu "num=1:j_exclusive=yes"
#BSUB -q gpu
#BSUB -m compgpu03
#BSUB -n 12

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load python/3.11.6
module load dorado/0.9.0

echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
export CUDA_VISIBLE_DEVICES

##reading in the variables:
pod5_files="/beevol/home/ouyangjo/raw_data/JiYoung_nanopore_data/FUET3Q0/YOU32634.08092025/YOU32634-s2_25AUG2025/250822_HeLa_Lipo_200uM/20250825_1331_1G_PBG54413_f72c934c/pod5/"
reference="/beevol/home/ouyangjo/genome_files/hg38_ENO2_rep1.fa"
output_folder="/beevol/home/ouyangjo/results/JiYoung/lipo_alone_FUET3Q0/hela_lipo_rep1_FUET3Q0/"
output_file="${output_folder%/}/hela_lipo_rep1_FUET3Q0.bam"
mkdir -p "$output_folder"

dorado basecaller \
  hac "$pod5_files" \
  --device cuda:0 \
  --estimate-poly-a \
  --reference "$reference" \
  --mm2-opts "-k 15 -w 10" \
  --min-qscore 9 \
  > "$output_file"
