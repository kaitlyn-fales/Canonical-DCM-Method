# Creation of VOI mask for Social Task DCM with HCP data

# Load package
library(RNifti)

# Initialize array of appropriate size
mask_right <- array(0, dim = c(91,109,91))

# Store coordinates for centers of VOIs
R_V5 <- c(24,32,39)
R_pSTS <- c(19,39,45)

# Add V5 into mask with a radius of 6mm (3 voxels) in each direction
mask_right[(R_V5[1]-3):(R_V5[1]+3),
           (R_V5[2]-3):(R_V5[2]+3),
           (R_V5[3]-3):(R_V5[3]+3)] <- 1

# Add pSTS into mask with a radius of 6mm (3 voxels) in each direction 
mask_right[(R_pSTS[1]-3):(R_pSTS[1]+3),
           (R_pSTS[2]-3):(R_pSTS[2]+3),
           (R_pSTS[3]-3):(R_pSTS[3]+3)] <- 2

# Initialize array of appropriate size
mask_left <- array(0, dim = c(91,109,91))

# Store coordinates for centers of VOIs
L_V5 <- c(68,27,39)
L_pSTS <- c(74,38,42)

# Add V5 into mask with a radius of 6mm (3 voxels) in each direction
mask_left[(L_V5[1]-3):(L_V5[1]+3),
          (L_V5[2]-3):(L_V5[2]+3),
          (L_V5[3]-3):(L_V5[3]+3)] <- 1

# Add pSTS into mask with a radius of 6mm (3 voxels) in each direction 
mask_left[(L_pSTS[1]-3):(L_pSTS[1]+3),
          (L_pSTS[2]-3):(L_pSTS[2]+3),
          (L_pSTS[3]-3):(L_pSTS[3]+3)] <- 2

# Load in sample data file (from HCP - might need to update paths)
#dat <- readNifti("100307_3T_tfMRI_SOCIAL_preproc/100307/MNINonLinear/Results/tfMRI_SOCIAL_RL/tfMRI_SOCIAL_RL.nii")

# Save mask and use data as template for header file
writeNifti(mask_right, "HCP_Social_Task/mask_R.nii", template = dat)
writeNifti(mask_left, "HCP_Social_Task/mask_L.nii", template = dat)
