# Automated Frequency Domain Decomposition (AFDD)
Automated Modal parameters identification from ambient vibrations measurement 


## Summary

The automated Frequency Domain Decomposition presented here is inspired by the Frequency Domain Decomposition (FDD) introduced by [1, 2]. The goal is to identify the mode shapes, eigenfrequencies and modal damping ratios from acceleration records obtained during structural health monitoring of civil engineering structures subjected to ambient noise. In this submission, an automated procedure is implemented in addition to the manual one proposed by [3]. For the automated procedure, I am using the peak picking function “pickpeaks” developed by [4] and available in [5], which was much more efficient than the Matlab function "findpeaks" for this purpose. I am, therefore, indebted to [3-5] for their previous works. The modal damping ratios are determined for each mode by using [6]. The acceleration data comes from a time-domain simulation of a clamped-free beam response to white noise excitation. The target modal properties from the beam come from [7].

## Content
The submission contains:

- The function AFDD

- an Matlab livescript file Documentation.mlx

- acceleration data beamData.m (4 Mb)

- The function pickpeaks.m [5]

Any comment, suggestion and question is welcome.

## References

[1] Brincker, R.; Zhang, L.; Andersen, P. (2001). "Modal identification of output-only systems using frequency domain decomposition". Smart Materials and Structures 10 (3): 441. doi:10.1088/0964-1726/10/3/303.

[2] Brincker, R., Zhang, L., & Andersen, P. (2000, February). Modal identification from ambient responses using frequency domain decomposition. In Proc. of the 18*‘International Modal Analysis Conference(IMAC), San Antonio, Texas.

[3] https://se.mathworks.com/matlabcentral/fileexchange/50988-frequency-domain-decomposition--fdd-

[4] Antoine Liutkus. Scale-Space Peak Picking. [Research Report] Inria Nancy - Grand Est (Villers-lès-Nancy, France). 2015. <hal-01103123v2>.
  
[5] https://se.mathworks.com/matlabcentral/fileexchange/42927-pickpeaks-v-select-display-

[6] https://se.mathworks.com/matlabcentral/fileexchange/55557-modal-parameters-identification-from-ambient-vibrations--sdof-

[7] https://se.mathworks.com/matlabcentral/fileexchange/52075-eigen-value-calculation-of-a-continuous-beam--transverse-vibrations-
