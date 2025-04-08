MRFMeanandVarmap.m calculates the mean and variance of MRF T1, MRF T2 values.
MRF_junctionmap_noCSFkernel.m generates MRF-based T1 and T2 z-score maps, The influence of CSF on MRF values were attenuated. MRF-based junction z-score maps were also generated based on MRF GM and WM tissue fractions. 
MRFextensionandthicknessmaps.m generates MRF-based extension z-score maps. The thickness maps were not generated because of low sensitivity.
generate_softmax.ipynb creates softmax outputs necessary to compare with MAP18 FCD probability outputs.
morphometry_own_cohort.m generates morphometric z-score maps based on MRF T1w image and our own healthy controls instead of using clinical T1w image with MAP18.
Analyze_nnUNet_results.ipynb calculates the performace metrics of outputs from Leave-one-out-cross validation of the model trained using nnU-Net.
nnunetsaliency.py generates saliency map though only of certain slices of certain layers
nnunet_numeric_saliency.py generates the averaged saliency of each input map at different layers across all models.
colorful_MRF_for_generating_figures.m generates all the figures used in the paper.