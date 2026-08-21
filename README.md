### Analysis workflow:
1. Run Cell Ranger v10 on raw FASTQs<br>
`cellranger_scripts/run_all.sh`
3. Determine downsampling fraction to deal with sequencing depth disparities<br>
`seqtk/get_downsample_fractions.R`
5. Run seqtk to downsample FASTQs appropriately<br>
`seqtk/seqtk.sh`
7. Re-run Cell Ranger v10 on downsampled FASTQs<br>
`cellranger_scripts/run_all_downsampled.sh`
9. Run Cell Bender<br>
`cellbender_scripts/run_all_cellbender.sh` 
11. Assemble objects<br>
`preprocessing/01_obj_creation.R`
