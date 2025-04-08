import torch
from batchgenerators.utilities.file_and_folder_operations import join
from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor
from nnunetv2.imageio.simpleitk_reader_writer import SimpleITKIO
import nibabel as nib
import numpy as np
from scipy import ndimage
import pandas as pd
import os
from pytorch_grad_cam import GradCAM, HiResCAM, ScoreCAM, GradCAMPlusPlus, AblationCAM, XGradCAM, LayerCAM, EigenCAM, FullGrad, GuidedBackpropReLUModel
from pytorch_grad_cam.utils.model_targets import ClassifierOutputTarget
from pytorch_grad_cam.utils.image import show_cam_on_image
from PIL import Image
import matplotlib.pyplot as plt
from PIL import Image as im
import matplotlib.cm as cm

subjIDs = [(studyID)]
slices = [70, 112 , 92, 108, 66, 97, 94, 127,
          88, 94, 34, 33, 56, 110, 122, 53,
          67, 111, 36, 53, 100, 123, 104, 117,
          129, 120, 132, 43, 107, 126, 68, 100,
          130, 84, 37, 81, 109, 117, 92, 73]
nnUNet_results = '/Volumes/eegrvw/Imaging/Multimodal/MRF/Peter/nnUNet_CV'
nnUNet_raw = '/Volumes/eegrvw/Imaging/Multimodal/MRF/Recon_MRF_3T/Patients/MRF_to_nnunet/nnUNet_raw'
predictor = nnUNetPredictor(tile_step_size=0.5,
                            use_gaussian=True,
                            use_mirroring=True,
                            perform_everything_on_device=True,
                            device=torch.device('cuda', 0),
                            verbose=False,
                            verbose_preprocessing=False,
                            allow_tqdm=True)
IDls = [(studyID)]

class SemanticSegmentationTarget:
    def __init__(self, category, mask):
        self.category = category
        self.mask = torch.from_numpy(mask)
        if torch.cuda.is_available():
            self.mask = self.mask.cuda()

    def __call__(self, model_output):
        return (model_output[:, :, :]).sum()

# initializes the network architecture, loads the checkpoint
df = pd.DataFrame(columns=['T1w', 'T1z', 'T2z', 'Juncz', 'Extenz', 'GM', 'WM', 'CSF'])
for stage in range(5,-1,-1):
    implists = []
    voxlists = []
    for fold, subjID in enumerate(subjIDs):
        if subjID in IDls:
            predictor.initialize_from_trained_model_folder(
                join(nnUNet_results, 'Dataset016_MRFz/nnUNetTrainer_250epochs__nnUNetPlans__3d_fullres'),
                use_folds=(fold,),
                checkpoint_name='checkpoint_best.pth')
            img, props = SimpleITKIO().read_images([join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0000.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0001.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0002.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0003.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0004.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0005.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0006.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0007.nii.gz'),
                                                    ])

            file = '/Volumes/eegrvw/Imaging/Multimodal/MRF/Peter/nnUNet_CV/Dataset016_MRFz/nnUNetTrainer_250epochs__nnUNetPlans__3d_fullres/fold_' + str(fold) + '/validation/' + subjID + '.nii.gz'
            output = nib.load(file).get_fdata()
            targets = [SemanticSegmentationTarget(1, output)]

            model = predictor.network.decoder.encoder.stages
            # model = predictor.network.encoder.stages

            target_layers = [predictor.network.decoder.encoder.stages[stage][0].convs[1]]
            cam = GradCAM(model=model, target_layers=target_layers)

            implist = []
            voxlist = []
            for chan in range(8):
                dup = np.zeros_like(img)
                dup[chan, :, :, :] = 1
                input_tensor = torch.from_numpy(img[:, :, :, :]*dup)
                input_tensor = input_tensor[None, :, :, :, :]
                grayscale_cam = cam(input_tensor=input_tensor, targets=targets)
                implist.append(grayscale_cam[grayscale_cam>0].sum())
                voxlist.append((grayscale_cam > 0).sum())
            implists.append(implist)
            voxlists.append(voxlist)

    implists = np.array(implists)
    voxlists = np.array(voxlists).astype(float)
    for i in range(implists.shape[0]):
        implists[i] = implists[i]/np.max(implists[i])
        voxlists[i] = voxlists[i]/np.max(voxlists[i])
    mimp = np.mean(implists, axis=0).transpose()
    mvox = np.mean(voxlists, axis=0).transpose()
    print(str(stage)+'imp: ' + str(mimp))
    print(str(stage) + 'vox: ' + str(mvox))
    ar1 = pd.DataFrame([mimp], columns=df.columns)
    ar2 = pd.DataFrame([mvox], columns=df.columns)
    df = pd.concat([df, ar1, ar2])
df.to_excel('/Users/irene/Desktop/featureimportance.xlsx', index=False)

df = pd.DataFrame(columns=['T1w', 'T1z', 'T2z', 'Juncz', 'Extenz', 'GM', 'WM', 'CSF'])
for stage in range(5,-1,-1):
    implists = []
    voxlists = []
    for fold, subjID in enumerate(subjIDs):
        if subjID in IDls:
            predictor.initialize_from_trained_model_folder(
                join(nnUNet_results, 'Dataset016_MRFz/nnUNetTrainer_250epochs__nnUNetPlans__3d_fullres'),
                use_folds=(fold,),
                checkpoint_name='checkpoint_best.pth')
            img, props = SimpleITKIO().read_images([join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0000.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0001.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0002.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0003.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0004.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0005.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0006.nii.gz'),
                                                    join(nnUNet_raw, 'Dataset016_MRFz/imagesTr/' + subjID +'_0007.nii.gz'),
                                                    ])

            file = '/Volumes/eegrvw/Imaging/Multimodal/MRF/Peter/nnUNet_CV/Dataset016_MRFz/nnUNetTrainer_250epochs__nnUNetPlans__3d_fullres/fold_' + str(fold) + '/validation/' + subjID + '.nii.gz'
            output = nib.load(file).get_fdata()
            targets = [SemanticSegmentationTarget(1, output)]

            model = predictor.network.encoder.stages
            # model = predictor.network.encoder.stages

            target_layers = [predictor.network.encoder.stages[stage][0].convs[1]]
            cam = GradCAM(model=model, target_layers=target_layers)

            implist = []
            voxlist = []
            for chan in range(8):
                dup = np.zeros_like(img)
                dup[chan, :, :, :] = 1
                input_tensor = torch.from_numpy(img[:, :, :, :]*dup)
                input_tensor = input_tensor[None, :, :, :, :]
                grayscale_cam = cam(input_tensor=input_tensor, targets=targets)
                implist.append(grayscale_cam[grayscale_cam>0].sum())
                voxlist.append((grayscale_cam > 0).sum())
            implists.append(implist)
            voxlists.append(voxlist)

    implists = np.array(implists)
    voxlists = np.array(voxlists).astype(float)
    for i in range(implists.shape[0]):
        implists[i] = implists[i]/np.max(implists[i])
        voxlists[i] = voxlists[i]/np.max(voxlists[i])
    mimp = np.mean(implists, axis=0).transpose()
    mvox = np.mean(voxlists, axis=0).transpose()
    print(str(stage)+'imp: ' + str(mimp))
    print(str(stage) + 'vox: ' + str(mvox))
    ar1 = pd.DataFrame([mimp], columns=df.columns)
    ar2 = pd.DataFrame([mvox], columns=df.columns)
    df = pd.concat([df, ar1, ar2])
df.to_excel('/Users/irene/Desktop/featureimportance_enco.xlsx', index=False)