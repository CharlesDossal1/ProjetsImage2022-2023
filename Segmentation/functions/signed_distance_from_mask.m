function phi= signed_distance_from_mask(mask)


Dm=fast_marching(mask);   %compute the distance to the region mask=1
Dp=fast_marching(1-mask); %compute the distance to the region mask=0
phi=Dp-Dm;