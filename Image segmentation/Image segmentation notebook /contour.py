
from roipoly import RoiPoly
import numpy as np
from matplotlib import pyplot as plt
from PIL import Image 

image =Image.open("/home/omezzine/Bureau/Image segmentation/coins.png").convert("L")
plt.imshow(image,cmap="gray")
initial_contour = RoiPoly(color='r')
contournp = initial_contour.get_mask(image).astype(int)
plt.figure()
plt.imshow(contournp)
plt.show()

np.save("Initial_contour",contournp)
initial_contour.display_roi()