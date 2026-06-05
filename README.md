# KineticFluxProfilingProject
## 	Simulation for k_app in the linear pathway
To get accurate estimation, we simulated the labeling fraction with different values of _r_ ($k_x/k_y$), and fit single exponential to obtain _k<sub>Y_app</sub>_. In this code, we compared the fitted _k<sub>app</sub>_ with the curve of the equation $k_{app}=k_X*k_Y/(k_X+k_Y)$, and the result shows that they are very similar.

![simulation for k_app in the linear pathway](https://github.com/user-attachments/assets/90dc7f4d-0de3-4b1f-a75f-ade4b9b3e735)

## 	Timepoint error simulation
To investigate the effect of time points on fitting error, we used this code for simulation. In this code, we sampled each time point under different values of k, calculated the fitting error arising from measurement error and biological error at that time point, and identified the optimal time point for each k value. These results can provide some guidance for our selection of time points.

<img width="1200" height="1000" alt="Figure_1" src="https://github.com/user-attachments/assets/ed10e45e-8657-441a-bd1b-05dacc88fb5d" />

<img width="1200" height="800" alt="Figure_2" src="https://github.com/user-attachments/assets/321ec38e-8735-40b6-8319-2dee01de8fd1" />

## 	Reversible reaction simulation for group criteria
