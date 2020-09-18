# Calibrated Intervention and Containment of the COVID-19 Pandemic

A MATLAB package for quantifying the COVID-19 pandemic through calibration with clinical data<sup>[1]</sup>

Liang Tian<sup>1,2,#</sup>, Xuefei Li<sup>1,3,#</sup>, Fei Qi<sup>1,3</sup>, Qian-Yuan Tang<sup>1,4</sup>, Viola Tang<sup>1,5</sup>, Jiang Liu<sup>1,6</sup>, Zhiyuan Li<sup>1,7</sup>, Xingye Cheng<sup>1,2</sup>, Xuanxuan Li<sup>1,8,9</sup>, Yingchen Shi<sup>1,8,9</sup>, Haiguang Liu<sup>1,9,10</sup>, Lei-Han Tang<sup>1,2,9,*</sup>
 
<sup>#</sup>Contributed equally      
<sup>*</sup>Correspondence: lhtang@hkbu.edu.hk
 
<sup>1</sup>COVID-19 Modelling Group, Hong Kong Baptist University, Kowloon, Hong Kong SAR, China, http://covid19group.hkbu.edu.hk  
<sup>2</sup>Department of Physics and Institute of Computational and Theoretical Studies, Hong Kong Baptist University, Kowloon, Hong Kong SAR, China  
<sup>3</sup>CAS Key Laboratory of Quantitative Engineering Biology, Shenzhen Institute of Synthetic Biology, Shenzhen Institutes of Advanced Technology, Shenzhen 518055, China  
<sup>4</sup>Center for Complex Systems Biology, Universal Biology Institute, University of Tokyo, 113-0033 Tokyo, Japan  
<sup>5</sup>Department of Information Systems, Business Statistics and Operations Management, Hong Kong University of Science and Technology, Hong Kong SAR, China
<sup>6</sup>Scarsdale, NY 10583, USA  
<sup>7</sup>Center for Quantitative Biology, Peking University, Haidian, Beijing 100871, China  
<sup></sup>Department of Engineering Physics, Tsinghua University, Haidian, Beijing 100084, China  
<sup>9</sup>Complex Systems Division, Beijing Computational Science Research Center, Haidian, Beijing 100193, China  
<sup>10</sup>Physics Department, Beijing Normal University, Haidian, Beijing 100875, China  

## Table of Contents
- [Download](#Download)
- [Files](#Files)
- [Installation](#Installation)
- [Outputs](#Outputs)
- [Reference](#Reference)

## Download   
<a href="https://github.com/hkbu-covid19group/Calib-covid19/blob/master/COVID19_tutorial_1.0.zip">COVID19_tutorial_1.0.zip</a>  (3.82 MB, version 1.0)  
Last update: Sep 17, 2020

## Files
The file COVID19_tutorial_1.0.zip contains:   
- Matlab scripts: Maximum likelihood estimation of symptom onset time distribution (Script_IP_MLE.m) and infectiousness profile (Script_Infec_MLE.m); Estimation of serial interval distribution through convolution (Script_SI_convolution.m); Calculation of basic reproduction number (Script_R0.m); Quantification of different intervention measurements, such as testing (Script_R0_reduction_testing.m), contact tracing (Script_R0_reduction_tracing.m), and mask-wearing (Script_R0_reduction_mask.m).
- Sample data: Symptom onset time of 347 confirmed cases (/Sampledata/IncubationPeriod.mat) from three independent studies: Bi et al., 183 cases<sup>[2]</sup>; Xia et al., 105 cases<sup>[3]</sup>; Men et al., 59 cases<sup>[4]</sup>. Exposure windows of 66 infector-infectee pairs (/Sampledata/ExposureWindow.mat)<sup>[5]</sup>.
- Sample figures: Main results.
- Readme: Instructions for this MATLAB package.

## Installation
- Extract the content of the enclosed COVID19_tutorial_1.0.zip file to a local directory.
- Run the corresponding MATLAB script files. Running times for the scripts are less than 1 minute.   

(The codes were written on MATLAB R2016b)

## Outputs
For details regarding the methods, please refer to our manuscript<sup>[1]</sup>.
- Symptom onset time distribution (Script_IP_MLE.m): Maximum likelihood estimation of the underlying symptom onset time distribution of the 347 cases data<sup>[2-4]</sup>. Estimated parameters mu, sigma, gamma and the corresponding likelihood are displayed in the command window.
![IP](https://github.com/hkbu-covid19group/Calib-covid19/blob/master/Samplefigures/IP.png)

- Infectiousness profile (Script_Infec_MLE.m): Maximum likelihood estimation of normalized infectiousness curve using the exposure window data of 66 infector-infectee pairs<sup>[5]</sup>. Estimated parameters alpha_A, alpha_B, theta_P, theta_S and the corresponding likelihood are displayed in the command window.
![Infec](https://github.com/hkbu-covid19group/Calib-covid19/blob/master/Samplefigures/Infec.png)

- Serial interval distribution inferred through convolution of symptom onset time distribution and infectiousness profile. The above estimated P_O(t) and P_I(t) are used as inputs.
![Infec](https://github.com/hkbu-covid19group/Calib-covid19/blob/master/Samplefigures/SI_conv.png)

- R0-lambda curve (Script_R0.m): Basic reproduction number R0 as a function of daily growth rate lambda. The above inferred symptom onset time distribution and infectiousness profile are used as inputs. Estimated R0 at lambda=0.3 is displayed in the command window.
![Infec](https://github.com/hkbu-covid19group/Calib-covid19/blob/master/Samplefigures/R0.png)

- Quantification of the effect of intervention measurements (Script_R0_reduction_testing.m, Script_R0_reduction_tracing.m): Calculation of reduced reproduction number Re under contact tracing (100% contacts are traced) and testing (no delay in testing report). The above inferred symptom onset time distribution and infectiousness profile are used as inputs.
![Infec](https://github.com/hkbu-covid19group/Calib-covid19/blob/master/Samplefigures/Intervention.png)

- Quantification of the effect of population-wide face mask wearing (Script_R0_reduction_mask.m): Left, reduction of transmission rate under mask wearing. Right, heatmap of reduced reproduction number Re as a function of percentage of mask-wearing population and percentage of traced contacts (at day 4 after infection). The above inferred symptom onset time distribution and infectiousness profile are used as inputs.
![Infec](https://github.com/hkbu-covid19group/Calib-covid19/blob/master/Samplefigures/mask-tracing.png)

## Reference
1.	Tian, L, et al. Calibrated intervention and containment of the COVID‐19 pandemic. arXiv 2003.07353 (2020).
2.	Xia, W. et al. Transmission of corona virus disease 2019 during the incubation period may lead to a quarantine loophole. medRxiv 2020.03.06.20031955 (2020). DOI:10.1101/2020.03.06.20031955
3.	Men K. et al. Estimate the incubation period of coronavirus 2019 (COVID-19). medRxiv 2020.02.24.20027474 (2020). DOI:10.1101/2020.02.24.20027474
4.	Bi, Q. et al. Epidemiology and Transmission of COVID-19 in Shenzhen China: Analysis of 391 cases and 1,286 of their close contacts. Lancet Infect. Dis. 20, 911-919, (2020)
5.	He, X. et al. Temporal dynamics in viral shedding and transmissibility of COVID-19. Nat. Med. 26, 672-675 (2020). 

