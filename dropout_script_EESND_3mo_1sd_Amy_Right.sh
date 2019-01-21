#!/bin/bash
# next line tells script to stop when there are undefined variables 
set -euo pipefail
IFS=$'\n\t'

STUDYDIR="/projects/sanchez_share/processed/EESND/"
FOLDERNAME="FAIRPRE_FSL_13_MONKEYINFANT3MOT1_DFM_NLIN"
ROI="Styner_to_F99_ROI_R.Amygdala_1.5mm_masked"
ROIDIR="/home/elyse.morin/testing/dropout/Amy_3mo_1sd_EESND_Right/"
GRAY_MATTER_ATLAS_MASK="/home/elyse.morin/testing/dropout/gray_matter_mask_in_F99_1.5mm_more_conservative.nii.gz"

#--- START PROCESSING LOOP
rm ${ROIDIR}subject_voxel_report
rm ${ROIDIR}${ROI}_Final_Union.nii.gz
full_ROI_voxels=$(cluster --in=${ROIDIR}${ROI}.nii.gz -t 1 | awk 'NR==2 {print $2}')
echo $ROI "includes" $full_ROI_voxels "voxels" >> ${ROIDIR}subject_voxel_report
echo "Subject" "threshold" "surviving_voxels" "coordinates" >> ${ROIDIR}subject_voxel_report
SUBNUM_INDEX=1
declare -i SUBNUM_INDEX

for SUBNUM in RFO14-3mo-T1rs-EM2 RKR14-3mo-T1rs-EM2 RWS14-3mo-T1rs-EM2 RNT14-3mo-T1rs-EM2 RST14-3mo-T1rs-EM2 RZU14-3mo-T1rs-EM2 RSA15-3mo-T1rs-EM4 RTA15-3mo-T1rs-EM2 RDF15-3mo-T1rs-EM2 RPH15-3mo-T1rs-EM2 RCK15-3mo-T1rs-EM2 RFN15-3mo-T1rs-EM2 RIC15-3mo-T1rs-EM2 RIO15-3mo-T1rs-EM2 RNC15-3mo-T1rs-EM2 ROH15-3mo-T1rs-EM2 RRF15-3mo-T1rs-EM2 RTM15-3mo-T1rs-EM2 RZK15-3mo-T1rs-EM2 RWS15-3mo-T1rs-EM2 RCT15-3mo-T1rs-EM2

do
cd ${ROIDIR}
mkdir -p ${SUBNUM}
cd ${SUBNUM}
cp ${STUDYDIR}${SUBNUM}/*/${FOLDERNAME}/REG*FUNC1/17*/xr3d_atl_res.4dfp.img ./ATLAS_EPI1.4dfp.img
cp ${STUDYDIR}${SUBNUM}/*/${FOLDERNAME}/REG*FUNC1/17*/xr3d_atl_res.4dfp.ifh ./ATLAS_EPI1.4dfp.ifh
caret_command -file-convert -vc ATLAS_EPI1.4dfp.ifh ATLAS_EPI1.nii.gz 
#Calculate Mean in gray matter mask, excluding voxels with all zero's
MEAN=`fslstats ATLAS_EPI1.nii.gz -k ${GRAY_MATTER_ATLAS_MASK} -M`
#Calculate Standard deviation in gray matter mask, excluding voxels with all zero's
STDEV=`fslstats ATLAS_EPI1.nii.gz -k ${GRAY_MATTER_ATLAS_MASK} -S`
#Calculate intensity threshold as mean minus (2 * Standard deviation) but with simple math
INTENSITY_THR=`echo "$MEAN-$STDEV" | bc`
# TURN INTO WHOLE NUMBER WITH NO DECIMAL
INTENSITY_THR=`echo "($INTENSITY_THR+0.5)/1" | bc`

fslmaths ATLAS_EPI1.nii.gz -mas ${ROIDIR}${ROI}.nii.gz ${ROI}_ATLAS_EPI1.nii.gz
##Threshold to the intensity threshold you created and binarize to create mask
fslmaths ${ROI}_ATLAS_EPI1.nii.gz -thr ${INTENSITY_THR} -bin ${ROI}_ATLAS_EPI1_bin.nii.gz
## Find number of voxels in subject-specific mask. Used cluster instead of fslstats -V because
# I was having trouble getting fslstats to report a reasonable number (maybe was inflating number of voxels by collapsing across time?)
# the NR==2 part isolsates the second line of the output. $2 isolates the second entry
voxels=$(cluster --in=${ROI}_ATLAS_EPI1_bin.nii.gz -t 1 | awk 'NR==2 {print $2}')

# Find coordinates that are within the range of the cut-off (50 above and 50 below) so you can go in and eyeball it to see if that cutoff makes sense
upper_thresh=`echo "$INTENSITY_THR+50" | bc`
lower_thresh=`echo "$INTENSITY_THR-50" | bc`

fslmaths ${ROI}_ATLAS_EPI1.nii.gz -thr $lower_thresh ${ROI}_ATLAS_EPI1_temp_lt.nii.gz
fslmaths ${ROI}_ATLAS_EPI1_temp_lt.nii.gz -uthr $upper_thresh ${ROI}_ATLAS_EPI1_temp_ut.nii.gz
fslmaths ${ROI}_ATLAS_EPI1_temp_ut.nii.gz -bin ${ROI}_ATLAS_EPI1_temp_ut_bin.nii.gz
cluster --in=${ROI}_ATLAS_EPI1_temp_ut_bin.nii.gz -t 1 > Cluster_report
output_Clust_report=$(wc -l "Cluster_report" | awk '{print $1}')

if [ $output_Clust_report == "1" ]; then
rm Cluster_report
upper_thresh=`echo "$INTENSITY_THR+100" | bc`
lower_thresh=`echo "$INTENSITY_THR-100" | bc`
fslmaths ${ROI}_ATLAS_EPI1.nii.gz -thr $lower_thresh ${ROI}_ATLAS_EPI1_temp_lt.nii.gz
fslmaths ${ROI}_ATLAS_EPI1_temp_lt.nii.gz -uthr $upper_thresh ${ROI}_ATLAS_EPI1_temp_ut.nii.gz
fslmaths ${ROI}_ATLAS_EPI1_temp_ut.nii.gz -bin ${ROI}_ATLAS_EPI1_temp_ut_bin.nii.gz
cluster --in=${ROI}_ATLAS_EPI1_temp_ut_bin.nii.gz -t 1 > Cluster_report
output_Clust_report=$(wc -l "Cluster_report" | awk '{print $1}')
fi

declare -i Output_Clust_report
if [ "$output_Clust_report" -gt "1" ]; then
coordinates=$(3dmaskdump -mask ${ROI}_ATLAS_EPI1_temp_ut.nii.gz -xyz ATLAS_EPI1.nii.gz | awk 'NR==1 {print "x=" $1 " y=" $2 " z="  $3}')
echo $SUBNUM $INTENSITY_THR $voxels "voxels" $coordinates >> ${ROIDIR}subject_voxel_report
else
echo $SUBNUM $INTENSITY_THR $voxels "no coordinates" >> ${ROIDIR}subject_voxel_report
fi

#if [ "$SUBNUM" = "RHd12_6dof" ]
if [ "$SUBNUM_INDEX" -eq 1 ]
then
cp ${ROI}_ATLAS_EPI1_bin.nii.gz  ${ROIDIR}${ROI}_Final_Union.nii.gz
SUBNUM_INDEX=$((SUBNUM_INDEX+1))
else
fslmaths ${ROI}_ATLAS_EPI1_bin.nii.gz -mul ${ROIDIR}${ROI}_Final_Union.nii.gz ${ROIDIR}${ROI}_Final_Union.nii.gz
fi
done
