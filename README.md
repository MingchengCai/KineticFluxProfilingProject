# KineticFluxProfilingProject
## 	Simulation for k_app in the linear pathway
To get accurate estimation, we simulated the labeling fraction with different values of _r_ ($k_x/k_y$), and fit single exponential to obtain _k<sub>Y_app</sub>_. In this code, we compared the fitted _k<sub>app</sub>_ with the curve of the equation $k_{app}=k_X*k_Y/(k_X+k_Y)$, and the result shows that they are very similar.

![simulation for k_app in the linear pathway](https://github.com/user-attachments/assets/90dc7f4d-0de3-4b1f-a75f-ade4b9b3e735)

## 	Timepoint error simulation
To investigate the effect of time points on fitting error, we used this code for simulation. In this code, we sampled each time point under different values of k, calculated the fitting error arising from measurement error and biological error at that time point, and identified the optimal time point for each k value. These results can provide some guidance for our selection of time points.

<img width="1200" height="1000" alt="Figure_1" src="https://github.com/user-attachments/assets/ed10e45e-8657-441a-bd1b-05dacc88fb5d" />

<img width="1200" height="800" alt="Figure_2" src="https://github.com/user-attachments/assets/321ec38e-8735-40b6-8319-2dee01de8fd1" />

## 	Reversible reaction simulation for group criteria
To verify the difference in fitting error between fitting two metabolites together versus fitting them separately in the presence of a reversible reaction, we conducted simulations using a simple reversible reaction structure. This allowed us to determine which fitting approach yields a smaller fitting error under different rate constants of the reversible reaction.

<img width="683" height="1138" alt="Figure_9" src="https://github.com/user-attachments/assets/34a866d4-ac93-4b48-98f7-e43bca07d33a" />

## 	Double EXP fitting
To fit metabolites that exhibit a significant delay and do not conform to a single exponential function (e.g., glutamine), we developed a code for double exponential fitting. In this code, we input the file k_value fitting.csv, in which every two metabolites are grouped as a pair, with the former metabolite considered as the direct upstream of the latter metabolite we aim to fit (e.g., glutamate and glutamine). The code performs double exponential fitting on the downstream metabolite by referencing the labeling rate of the upstream metabolite.
