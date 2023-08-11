# Automated Frequency Domain Decomposition (AFDD)
Automated Modal parameters identification from ambient vibrations measurement 

[![View Automated Frequency Domain Decomposition (AFDD) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/57153-automated-frequency-domain-decomposition-afdd)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4277622.svg)](https://doi.org/10.5281/zenodo.4277622)
<a href="https://www.buymeacoffee.com/echeynet" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/orange_img.png" alt="Buy Me A Coffee" style="height: 25px !important;width: 120px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a>
## Summary
The automated Frequency Domain Decomposition presented was applied in [1]. It inspired by the Frequency Domain Decomposition (FDD) introduced by [2, 3]. The goal is to identify the mode shapes, eigenfrequencies and modal damping ratios from acceleration records obtained during structural health monitoring of civil engineering structures subjected to ambient noise. In this submission, an automated procedure is implemented in addition to the manual one proposed by [4]. For the automated procedure, I am using the peak picking function “pickpeaks” developed by [5] and available in [6], which was much more efficient than the Matlab function "findpeaks" for this purpose. I am, therefore, indebted to [4-6] for their previous works. The modal damping ratios are determined for each mode by using [7]. The acceleration data comes from a time-domain simulation of a clamped-free beam response to white noise excitation. The target modal properties from the beam come from [8].

## Content
The submission contains:

- The function AFDD
- A Matlab livescript file Documentation.mlx
- Acceleration data beamData.m (4 Mb)
- The function pickpeaks.m [6]

Any comment, suggestion and question is welcome.

## References

[1] Cheynet, E., Jakobsen, J. B., & Snæbjörnsson, J. (2017). Damping estimation of large wind-sensitive structures. Procedia engineering, 199, 2047-2053.

[2] Brincker, R.; Zhang, L.; Andersen, P. (2001). "Modal identification of output-only systems using frequency domain decomposition". Smart Materials and Structures 10 (3): 441. doi:10.1088/0964-1726/10/3/303.

[3] Brincker, R., Zhang, L., & Andersen, P. (2000, February). Modal identification from ambient responses using frequency domain decomposition. In Proc. of the 18*‘International Modal Analysis Conference(IMAC), San Antonio, Texas.

[4] https://se.mathworks.com/matlabcentral/fileexchange/50988-frequency-domain-decomposition--fdd-

[5] Antoine Liutkus. Scale-Space Peak Picking. [Research Report] Inria Nancy - Grand Est (Villers-lès-Nancy, France). 2015. .

[6] https://se.mathworks.com/matlabcentral/fileexchange/42927-pickpeaks-v-select-display-

[7] https://se.mathworks.com/matlabcentral/fileexchange/55557-modal-parameters-identification-from-ambient-vibrations--sdof-

[8] https://se.mathworks.com/matlabcentral/fileexchange/52075-eigen-value-calculation-of-a-continuous-beam--transverse-vibrations-
