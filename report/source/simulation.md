# Point Simulation
To conduct a preliminary study and test the camera pose optimization and planning algorithm, the problem is recreated in a simplified manner within a synthetic, controlled, and reproducible environment to create, using a procedural method, the largest dataset of configurations possible. In this phase, the change detection component of the pipeline is assumed to be given and perfectly working; therefore, only the components using the point cloud are reproduced. The photometric difference is momentarily ignored and reduced to a proxy.

## Synthetic Point Cloud
A real depth sensor (RGB-D) or a structure from motion (SfM) outputs only the reconstruction of the scene's surfaces; it cannot provide the points contained inside the volume of the present objects. This has direct consequences for reconstruction properties, such as occlusions, incident angles, and the boundaries of known and unknown zones. To simulate this, it is necessary to first generate the surface of a random volume and then sample its points.

Using a level surface, it is possible to define a function $g(x,y,z)$ associating a number each point of 3D space:
- $g(x,y,z)<0$, the point is inside the solid volume;
- $g(x,y,z)<0$, the point is on the surface of the volume;
- $g(x,y,z)>0$, the is point outside the solid volume.


We define a random surface using placeholder formulation 



`\ref{...}`{=tex}