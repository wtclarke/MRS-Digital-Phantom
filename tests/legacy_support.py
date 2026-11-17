"""Small controlled anatomy for legacy regression tests."""
import numpy as np
import torch
import torchio as tio


def anatomy():
    affine=np.diag([2.,2.,2.,1.])
    labels=torch.tensor([0,1,1,2,2,2,3,3],dtype=torch.int16).reshape(1,2,2,2)
    return tio.Subject(labels=tio.LabelMap(tensor=labels,affine=affine),
        B0=tio.ScalarImage(tensor=torch.full((1,2,2,2),.5),affine=affine),
        brain_mask=tio.LabelMap(tensor=torch.ones(1,2,2,2),affine=affine))

