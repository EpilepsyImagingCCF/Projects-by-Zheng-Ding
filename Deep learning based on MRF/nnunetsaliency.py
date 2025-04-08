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

def apply_colormap(data, cmap_name='bwr'):
    cmap = cm.get_cmap(cmap_name)
    normed_data = (data - np.min(data)) / (np.max(data) - np.min(data))  # Normalize data to 0-1
    colormapped_data = cmap(normed_data)  # Apply colormap
    return (colormapped_data[:, :, :3]*255).astype(np.uint8)  # Convert to 8-bit RGB

class SemanticSegmentationTarget:
    def __init__(self, category, mask):
        self.category = category
        self.mask = torch.from_numpy(mask)
        if torch.cuda.is_available():
            self.mask = self.mask.cuda()

    def __call__(self, model_output):
        print((model_output[:, :, :]).sum())
        return (model_output[:, :, :]).sum()

# initializes the network architecture, loads the checkpoint

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

        # model = predictor.network.decoder
        # # target_layers = [predictor.network.encoder.stages[0][0].convs[0]]
        # target_layers = [predictor.network.decoder.stages[0]]
        # cam = GradCAM(model=model, target_layers=target_layers)
        # input_tensor = torch.from_numpy(img)
        # grayscale_cam = cam(input_tensor=input_tensor, targets=targets)
        # print(grayscale_cam.shape)

        model = predictor.network.decoder.encoder.stages
        model = predictor.network.encoder.stages

        # data = im.fromarray((img[0, 112, :, :]*255).squeeze().astype(np.uint8))
        # data.save("/Users/irene/Desktop/figure1.png")
        #
        # data = im.fromarray((grayscale_cam[0, 112, :, :]*255).squeeze().astype(np.uint8))
        # data.save("/Users/irene/Desktop/figure2.png")

        # Image.fromarray(cam_image)
        # fig = plt.imshow((grayscale_cam[0, 100, :, :].squeeze()).astype(np.uint8))
        # fig.savefig('/Users/irene/Desktop/figure.png')
        # print(grayscale_cam[0, 100, :, :].squeeze())
        # plt.imshow(grayscale_cam[100,:,:])

        # Function to apply a colormap to a numpy array

        zslice = slices[fold]
        for stage in range(4,5):
            target_layers = [predictor.network.encoder.stages[stage][0].convs[1]]
            cam = GradCAM(model=model, target_layers=target_layers)

            fig, axs = plt.subplots(2, 4, figsize=(40, 10))
            for chan in range(8):
                dup = np.zeros_like(img)
                dup[chan, :, :, :] = 1
                input_tensor = torch.from_numpy(img[:, :, :, :]*dup)
                input_tensor = input_tensor[None, :, :, :, :]
                grayscale_cam = cam(input_tensor=input_tensor, targets=targets)

                grayscale_image = im.fromarray((img[chan, zslice, :, :]*255).squeeze().astype(np.uint8))
                if (chan>0) and (chan<5):
                    probability_image = (grayscale_cam[0, zslice, :, :] * 25).squeeze().astype(np.uint8)
                else:
                    probability_image = (grayscale_cam[0, zslice, :, :] * 255).squeeze().astype(np.uint8)
                colormapped_prob_image = apply_colormap(probability_image)
                colormapped_image_pil = Image.fromarray(colormapped_prob_image)

                overlay_image = Image.blend(grayscale_image.convert('RGB'), colormapped_image_pil, alpha=0.3)
                overlay_image = overlay_image.rotate(180)
                # overlay_image.save('/Users/irene/Desktop/overlay_image.png')
                axs[int(chan/4), chan%4].imshow(overlay_image)
                axs[int(chan/4), chan%4].axis('off')

            # input_tensor = torch.from_numpy(img[:, :, :, :])
            # input_tensor = input_tensor[None, :, :, :, :]
            # grayscale_cam = cam(input_tensor=input_tensor, targets=targets)
            #
            # grayscale_image = img
            # cam_image = show_cam_on_image(grayscale_image, grayscale_cam, use_rgb=True)
            # cam_image = camimage[0,zslice,:,:]
            # cam_image.save('/Users/irene/Desktop/layer' + str(stage) + '.png')

            # uncomment this do do overall instead of
            # input_tensor = torch.from_numpy(img[:, :, :, :])
            # input_tensor = input_tensor[None, :, :, :, :]
            # grayscale_cam = cam(input_tensor=input_tensor, targets=targets)
            # grayscale_image = im.fromarray((img[0, zslice, :, :]*255).squeeze().astype(np.uint8))
            # probability_image = (grayscale_cam[0, zslice, :, :]*255).squeeze().astype(np.uint8)
            # colormapped_prob_image = apply_colormap(probability_image)
            # colormapped_image_pil = Image.fromarray(colormapped_prob_image)
            # overlay_image = Image.blend(grayscale_image.convert('RGB'), colormapped_image_pil, alpha=0.3)
            # overlay_image = overlay_image.rotate(180)
            # overlay_image.save('/Users/irene/Desktop/layercombined' + str(stage) + '.png')

            fig.tight_layout()
            fig.canvas.draw()
            image_from_plot = Image.frombytes('RGB', fig.canvas.get_width_height(), fig.canvas.tostring_rgb())
            image_from_plot.save('/Users/irene/Desktop/layer' + str(stage) + '.png')


        # image1 = im.fromarray((img[0, 112, :, :]*255).squeeze().astype(np.uint8))
        # image2 = im.fromarray((img[0, 112, :, :]*255).squeeze().astype(np.uint8))
        # image3 = im.fromarray((img[0, 112, :, :]*255).squeeze().astype(np.uint8))
        # image4 = im.fromarray((img[0, 112, :, :]*255).squeeze().astype(np.uint8))
        #
        # # Create a figure to display the images
        #
        #
        # # Display the images in the figure

        # axs[0, 1].imshow(image2, cmap='gray')
        # axs[0, 1].axis('off')
        # axs[1, 0].imshow(image3, cmap='gray')
        # axs[1, 0].axis('off')
        # axs[1, 1].imshow(image4, cmap='gray')
        # axs[1, 1].axis('off')
        #
        # # # Save the figure to a PIL Image
        # # fig.canvas.draw()
        # # image_from_plot = Image.frombytes('L', fig.canvas.get_width_height(), fig.canvas.tostring_grayscale())
        # #
        # # # Save the image
        # # image_from_plot.save('/Users/irene/Desktop/multi_image_plot.png')

        #
        # # Save the plot to a buffer
        # fig.canvas.draw()
        # buf = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
        # buf = buf.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        #
        # # Convert buffer to PIL Image
        # image = Image.fromarray(buf)
        #
        # # Save image
        # image.save('/Users/irene/Desktop/multi_grayscale_plot.png')