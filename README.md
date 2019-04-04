# FNN (false nearest neighbors) embedding
False nearest neighbors (FNN) method is proposed by Kennel(Kennel, 1992) to find minimal embedding dimension for time series dynamic systems. When applied to time series neural data, the minimal embedding dimension can reflect the complexity of neural activities and may provide a way to indicate disease states. Here we tested FNN method on two syhthetic systems. One is Henon map which is a two parameters dynamical system that exhibit chaotic behavior, the other one is synthetic electrode signals of neurons(Smith, 2019). 

1. Henon map: run 'FNN_HenonMap.m'. The minimal embedding dimension is two, which is consistent with the two parameters. When the number of time point is small (say 1000), the FNNP would increase at large dimension, increasing the number of time points can reduce this effect.

2. Synthetic eletrode signals: run 'FNN_Electrode.m'. The minimal embedding dimension of signals with 2kHz sampling rate is higher than the signals with 10kHz sampling rate, this may be explained by the fact that downsampling of singals reduces the information and increases the complexity of the system. 

## Reference
M. B. Kennel, R. Brown, and H. D. I. Abarbanel, Determining embedding dimension for phase-space reconstruction using a geometrical construction, Phys. Rev. A 45, 3403 (1992).

Smith, Leslie & Mtetwa, Nhamoinesu. (2019). Manual for the noisy spike generator matlab software. 
