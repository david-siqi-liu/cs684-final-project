# CS684 Final Project
## Goals

Train on labeled images in sketch_train.txt (#49,115) /real_train.txt (#122,563) /quickdraw_train.txt (#120,750) /infograph_train.txt (#37,087)

Test on unlabeled images in clipart_test.txt (#14,814).

Modifications:

1. 10 labels

## Models

Naive, Pre-trained model, transfer learning:

- ResNet-34
- No transformations (except for resize, normalization and to tensor)

Transformations:

- Test: none
- Train (one at a time)
  - Center-crop
  - Random-crop
  - Horizontal flips
  - Random rotation
  - Horizontal shifts
  - Vertical shifts
  - Color Jitter
  - RandomAffineTransformation

Augmentations (adding channels)

- Test: none
- Train 
  - Greyscale layer
  - Horizontal gradient (canny edge detector)
  - Vertical gradient

Conditional auto-encoder?

Add:

- Atrous
- Skip

