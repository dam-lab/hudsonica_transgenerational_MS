# hudsonica_transgenerational_MS
This repository is a collection of code and data for the manuscript titled: Limited copepod adaptation to combined warming and acidification reveals cost of producing adaptive phenotypes

This is a README file for the repository linked to the manuscript titled “Limited copepod adaptation to combined warming and acidification reveals cost of producing adaptive phenotypes”. There are two directories within the repository. This document will outline the contents of each directory and the intended use for each.

I.	Code
The first directory named “Code” contains all R scripts used to analyze and visualize the data presented in the manuscript. There are six (6) total scripts used, which are outlined below:

1.	A_hudsonica_physical_data.R - 
This script analyzes the physico-chemical data collected during the experiment from June 2019-August 2020. The script reads available chronological data (see “Data” section below) of temperature, pH, and pCO2, and other chemical parameters calculated from alkalinity titration experiments. This script evaluates whether physico-chemical data is similar across treatments.

2.	Body_size_analysis_MS_hudsonica.R - 
This script reads in Body size data and tests for changes in body size across generations for each treatment. This script also calculates Growth rate data based on Development time data.

3.	EPR_complete_MS_hudsonica.R - 
This script evaluates egg production rate (EPR) and hatching success (HS) data. It tests for changes in trait values across generations for each treatment and tests for effects of temperature/pH on EPR and HS. It also shows the scripts used to reproduce the figures presented in the manuscript.

4.	Fitness_complete_MS_hudsonica.R - 
This script calculates fitness (λ) data based on available survival, EPR, HS, sex ratio, and development time data. This script tests for changes in λ across generations for each treatment and evaluates the contributions of changing life-history traits to relative fitness. It also has scripts used to reproduce the figures presented in the manuscript.

5.	Reciprocal_transplant_data_MS_hudsonica.R - 
This script evaluates all survival, EPR, HS, sex ratio, development time, and λ data for experiments completed during the reciprocal transplant portion of the experiment. It tests for lineage X environment interactions and has scripts used to reproduce the figures in the manuscript.

6.	SurvivalData_complete_MS_hudsonica.R -
This script evaluates survival data collected during the experiment. It tests for changes in survival across generations for each treatment and includes scripts for reproducing figures presented in the manuscript.

7.	fst_for_Ahudsonica_MS.R
This script calculates Fst values for each lineage at each generation relative to the AM lineage for the same generation. Uses the “filtered_variants_Ahudsonica_MS.txt” file and the variants sync file output by “to_sync_Ahudsonica_MS.py”. Produces the “fst_Ahudsonica_MS.txt” file for later visualization.

8.	to_sync_Ahudsonica_MS.py
This script converts the filtered_variants_Ahudsonica_MS.txt file to sync format. 


II.	Data
The second directory named “Data” contains all data files needed to be analyzed including some output files generated from the scripts. There are a total of eleven (11) files. For clarity, files needed for analysis will be described first. Output data files will be described after.

A.	Input data files

1.	Physical_data_MS.xlsx – this file is a Microsoft Excel Workbook that contains fourteen (14) tabs of summary physico-chemical data. Tabs in the workbook are:
a)	Alkalinity_data_total_leuker2000 – a summary of all the calculated carbonate chemistry metrics from CO2SYS based on measured temperature, pH, salinity, and total alkalinity. This is the same table as Alkalinity_data_total_leuker2000.txt.
b)	Physical_data_MS – a chronological table of pH, Temperature, and input pCO2. **Note: the input pCO2 is different than the actual, measured value of pCO2 calculated during titrations. This is the same table as Physical_data_MS.txt.
c)	Temp mean – a summary table of the measured temperatures with standard deviations (sd), number of observations (n), standard error (se), and 95% confidence intervals (lower.ci & upper.ci).
d)	pCO2 mean – a summary table of the measured input pCO2 values with standard deviations (sd), number of observations (n), standard error (se), and 95% confidence intervals (lower.ci & upper.ci).
e)	 pH mean – a summary table of the measured pH values with standard deviations (sd), number of observations (n), standard error (se), and 95% confidence intervals (lower.ci & upper.ci).
f)	Temp contrasts – a summary table of contrasts of temperature measurements between treatments (including reciprocal transplants). Contrasts are noted in the first column named ‘contrast’ where the two treatments evaluated are separated by a hyphen (‘-’).
g)	pH contrasts – a summary table of contrasts of pH measurements between treatments (including reciprocal transplants). Contrasts are noted in the first column named ‘contrast’ where the two treatments evaluated are separated by a hyphen (‘-’).
h)	pCO2 contrasts – a summary table of contrasts of pCO2 measurements between treatments. Contrasts are noted in the first column named ‘contrast’ where the two treatments evaluated are separated by a hyphen (‘-’).
i)	TA mean – a summary table of the measured total alkalinity (TA) values with standard deviations (sd), number of observations (n), standard error (se), and 95% confidence intervals (lower.ci & upper.ci).
j)	mean salinity – a summary table of the mean measured salinity values.
k)	omega CA mean – a summary table of the measured omega calcium (ΩCA) values with standard deviations (sd), number of observations (n.count), standard error (se), and 95% confidence intervals (lower.ci & upper.ci).
l)	omega AR mean – a summary table of the measured omega aragonite (ΩAR) values with standard deviations (sd), number of observations (n.count), standard error (se), and 95% confidence intervals (lower.ci & upper.ci).
m)	fCO2 mean – a summary table of the measured fugacity (fCO2) values with standard deviations (sd), number of observations (n.count), standard error (se), and 95% confidence intervals (lower.ci & upper.ci).
n)	DIC mean – a summary table of the dissolved inorganic carbon values with standard deviations (sd), number of observations (n.count), standard error (se), and 95% confidence intervals (lower.ci & upper.ci).

2.	Alkalinity_data_total_leuker2000.txt – this file is a summary of all the calculated carbonate chemistry metrics from CO2SYS based on measured temperature, pH, salinity, and total alkalinity. This file is read into the script “A_hudsonica_physical_data.R”.

3.	Physical_data_MS.txt – this file is a chronological table of pH, Temperature, and input pCO2. **Note: the input pCO2 is different than the actual, measured value of pCO2 calculated during titrations. This file is read into the script “A_hudsonica_physical_data.R”.

4.	Survival_data_total.txt – This file is a tab-delimited text file of raw data used to calculate survival and development rates. There are thirteen (13) columns of data:
a)	Generation – the generation where the data was collected.
b)	Treatment – the treatment affiliated with survival data collected (1 = AM, 2 = OA, 3 = OW, 4 = OWA).
c)	Temp – the temperature (°C) of the respective treatment.
d)	pH – the pH of the respective treatment.
e)	Rep – the culture replicate that the data was collected from.
f)	Beak – the beaker used during the experiment.
g)	time – the day evaluated for the survival experiment.
h)	nx – the number of live individuals on any day.
i)	lx - the proportion of surviving individuals.
j)	Ndev – the number of copepodites observed on any day.
k)	Cdev – the number of adults observed on any day.
l)	F.Ratio – the ratio of females observed on any day.
m)	M.Ratio – the ratio of females observed on any day.

5.	EPR_HS_data_total_raw.txt – This file is a tab-delimited text file of raw data used to calculate egg production rate and hatching success. There are a total of seventeen (17) columns of data:
a)	Generation – the generation where the data was collected.
b)	Number – a replicate number assigned to an individual mate pair evaluated during egg production experiments.
c)	Treatment – the treatment affiliated with survival data collected.
d)	Temp – the temperature (°C) of the respective treatment.
e)	pH – the pH of the respective treatment.
f)	Naup1 – the number of nauplii observed after the first 48 h egg laying period.
g)	Un1 – the number of unhatched eggs observed after the first 48 h egg laying period.
h)	Sum1 – the sum of Naup1 and Un1.
i)	EPR1 – the egg production rate over the second 48 h.
j)	HF1 – the hatching success of eggs laid over the second 48 h.
k)	Naup2 – the number of nauplii observed after the second 48 h egg laying period.
l)	Un2 – the number of unhatched eggs observed after the second 48 h egg laying period.
m)	Sum2 – the sum of Naup1 and Un1.
n)	EPR2 – the egg production rate over the second 48 h.
o)	HF2 – the hatching success of eggs laid over the second 48 h.
p)	EPRtot – the overall rate of egg production over the entire 96 h.
q)	Hftot – the overall rate of hatching success.

6.	Feeding_cost_epr.txt – this file is a tab-delimited text file of egg production and hatching success data collected during the F11 reciprocal transplant food-limitation experiment. There are a total of fourteen (14) columns:
a)	Treatment – the treatment affiliated with survival data collected (1 = AM, 2 = OA, 3 = OW, 4 = OWA).
b)	Food – the food concentration (μg C/L) for each individual mate pair.
c)	Bin – the culture replicate that the data was collected from.
d)	Hatch1 – the number of hatched nauplii observed after the first 24 h egg laying period.
e)	Hatch2 – the number of hatched nauplii observed after the second 24 h egg laying period.
f)	Hatch3 – the number of hatched nauplii observed after the third 24 h egg laying period.
g)	Unhatch1 – the number of unhatched eggs observed after the first 24 h egg laying period.
h)	Unhatch2 – the number of unhatched eggs observed after the second 24 h egg laying period.
i)	Unhatch3 – the number of unhatched eggs observed after the third 24 h egg laying period.
j)	EPR – the overall rate of egg production over the entire 72 h egg laying period.
k)	HF1 – the hatching frequency after the first 24 h egg laying period.
l)	HF2 – the hatching frequency after the second 24 h egg laying period.
m)	HF2 – the hatching frequency after the third 24 h egg laying period.
n)	HF – the overall hatching frequency over the entire 72 h egg laying period.

7.	SurvDataTransplant.txt – this file is a tab-delimited text file of survival data collected during the F11 reciprocal transplant food-limitation experiment. There are a total of thirteen (13) columns:
a)	Treatment – the treatment affiliated with survival data collected (1 = AM, 4 = OWA, 5 = OWA->AM, 6 = AM->OWA).
b)	Temp – the temperature (°C) of the respective treatment.
c)	pH – the pH of the respective treatment.
d)	Rep – the culture replicate that the data was collected from.
e)	Beak – the beaker used during the experiment.
f)	time – the day evaluated for the survival experiment.
g)	nx – the number of live individuals on any day.
h)	lx - the proportion of surviving individuals.
i)	Ndev – the number of copepodites observed on any day.
j)	Cdev – the number of adults observed on any day.
k)	F.Ratio – the ratio of females observed on any day.
l)	M.Ratio – the ratio of females observed on any day.

8.	Body_size_data_MS_complete.txt – this file is a tab-delimited text file of body size data collected for the four original treatments at generations 0, 2, 4, and 11. The F11 data only includes the AM and OWA treatments for the food-limited experiment. There are six (6) columns of data:
a)	Generation – the generation where the data was collected.
b)	Treatment – the treatment affiliated with survival data collected (1 = AM, 2 = OA, 3 = OW, 4 = OWA).
c)	Replicate – the culture replicate that the data was collected from.
d)	Stage – the life stage that the individual was collected and measured at (C1 = C1, C6F = C6 Female, C6M = C6 Male.
e)	Number – the individual within a replicate for each treatment.
f)	Length – the measured length (mm).

9.	filtered_variants_Ahudsonica_MS.txt – this is a tab-delimited text file for evaluating Fst values. There are 48 columns. Details on column information can be found at https://varscan.sourceforge.net/using-varscan.html. 

B.	Output Data files

1.	EPR_HF_data_total_w_11.txt – this file is a tab-delimited text file of data used to calculate egg production rate (EPR) and hatching success (HS) for all generations of the original four treatments (AM, OA, OW, and OWA). It is summarized from the EPR_HS_data_total_raw.txt and Feeding_cost_epr.txt files. There are nine (9) columns of data:
a)	Generation – the generation where the data was collected.
b)	Number – a replicate number assigned to an individual mate pair evaluated during egg production experiments.
c)	Treatment – the treatment affiliated with survival data collected (1 = AM, 2 = OA, 3 = OW, 4 = OWA).
d)	Temp – the temperature (°C) of the respective treatment.
e)	pH – the pH of the respective treatment.
f)	EPRtot – the calculated egg production rate.
g)	Hftot – the calculated hatching success rate.
h)	Generation.c – the generation where the data was collected.
i)	Rep – the culture replicate that the data was collected from.

2.	Development_time_w_11.txt – this file is a summarized tab-delimited text file of mean development time data calculated from the raw survival data. The F11 data only includes the AM and OWA treatments for the food-limited experiment. There are eleven (11) columns of data:
a)	Treatment – the treatment affiliated with survival data collected (1 = AM, 2 = OA, 3 = OW, 4 = OWA).
b)	Rep – the culture replicate that the data was collected from.
c)	Beak – the beaker used during the experiment.
d)	Generation – the generation where the data was collected.
e)	mean – the mean calculated development time.
f)	sd – the standard deviation of the calculated development time.
g)	n.count – the number of individuals used to calculate development time.
h)	se – the standard error of the calculated development time.
i)	lower.ci – the lower bound of the 95% confidence interval.
j)	upper.ci – the upper bound of the 95% confidence interval.
k)	Stage – the life stage reflecting the calculated development time.

3.	GR_data_MS_w_11.txt – this file is a summary of calculated somatic growth rates from measured body sizes and calculated development times. The F11 data only includes the AM and OWA treatments for the food-limited experiment. There are fifteen (15) columns of data:
a)	Generation – the generation where the data was collected.
b)	Treatment – the treatment affiliated with survival data collected (1 = AM, 2 = OA, 3 = OW, 4 = OWA).
c)	Rep – the culture replicate that the data was collected from.
d)	Stage.x – the stage of the mature life stage used to calculate the somatic growth rate.
e)	Number.x – the individual within a replicate for each treatment of the mature life stage.
f)	Length.x – the measured length of the mature individual.
g)	Weight.x – the calculated weight of the mature individual.
h)	Stage.y – the stage of the C1 individual used to pair with the mature individual for change in body weight.
i)	Number.y – the individual within a replicate for each treatment of the C1 life stage.
j)	Weight.y – the calculated weight of the C1 individual.
k)	weight.diff – the change in weight between the C1 and C6F/C6M stage.
l)	Beak – the beaker used during the experiment.
m)	Dev.diff – the difference in development time between the mature life stage and the C1 life stage.
n)	Growth.Rate – the calculated somatic growth rate.

4.	lambda_results_devtime_surv_epr_hf_sex_w_f11.txt – This is a tab-delimited text file of the calculated fitness (λ) data derived from life-history trait data. The F11 data only includes the AM and OWA treatments for the food-limited experiment. There are nine (9) columns of data:
a)	Generation – the generation where the data was collected.
b)	Treatment – the treatment affiliated with survival data collected (1 = AM, 2 = OA, 3 = OW, 4 = OWA).
c)	Rep – the culture replicate that the data was collected from.
d)	lambda – the calculated fitness (λ) value.
e)	surv – the survival value used to calculate λ.
f)	epr – the EPR value (eggs per female per day) used to calculate λ.
g)	hf – the hatching success value used to calculate λ.
h)	sex – the female sex ratio used to calculate λ.
i)	dev.time – the development time (days) to adulthood needed to calculate λ.

5.	lambda_results_cost_f11_surv_epr_hf_devtime.txt – This is a tab-delimited text file of the calculated fitness (λ) data derived from life-history trait data during the F11 reciprocal transplant and food-limitation experiments. There are eight (8) columns of data:
a)	Treatment – the treatment affiliated with survival data collected (1 = AM, 4 = OWA, 5 = OWA->AM, 6 = AM->OWA).
b)	Food – the food concentration (μg C/L) for each individual mate pair.
c)	Rep – the culture replicate that the data was collected from.
d)	surv – the survival value used to calculate λ.
e)	epr – the EPR value (eggs per female per day) used to calculate λ.
f)	hf – the hatching success value used to calculate λ.
g)	dev.time – the development time (days) to adulthood needed to calculate λ.
h)	lambda – the calculated fitness (λ) value.

6.	fst_Ahudsonica_MS.txt – this is a tab-delimited file of Fst values calculated from variant sites. There are 14 columns of data:
a)	Samp1 – the first sample being compared.
b)	Treatment_1 – the treatment of the first sample being compared.
c)	Generation_1 – the generation of the first sample being compared.
d)	Rep_1 – the replicate of the first sample being compared.
e)	Samp2 – the second sample being compared.
f)	Treatment_2 – the treatment of the second sample being compared.
g)	Generation_2 – the generation of the second sample being compared.
h)	Rep_2 – the replicate of the second sample being compared.
i)	fst – the calculated Fst value.
j)	trt_gen1 – the treatment and generation combined for the first sample being combined.
k)	trt_gen2 – the treatment and generation combined for the second sample being combined.
l)	group – the two treatments being compared.
m)	group_gen – the two treatments being compared with the relevant generation.
n)	generation_mod – the generation for the two treatments being compared.


III.	Statistics

A.	Summary files – all files in this directory are MS Excel Workbooks with multiple tabs that consist of model summaries and ANOVA tables. All tabs are labeled with the respective model type for a particular trait. 

1.	BodySize_GrowthRate_stats_MS_updated.xls – this is a spreadsheet with summarized statistics for Body Size and Growth Rate across the transgenerational experiment.
2.	Dev_time_stats_MS_updated.xls – this is a spreadsheet with summarized statistics for Development time across the transgenerational experiment.
3.	EPR_HS_stats_MS.xls – this is a spreadsheet with summarized statistics for EPR and HS across the transgenerational experiment.
4.	EPR_RT_statistics_MS.xls – this is a spreadsheet with summarized statistics for EPR in the reciprocal transplant part of the experiment.
5.	Lambda_RT_statistics_MS.xls – this is a spreadsheet with summarized statistics for calculated fitness in the reciprocal transplant part of the experiment.
6.	Lambda_stats_MS_updated.xls – this is a spreadsheet with summarized statistics for calculated fitness across the transgenerational experiment.
7.	Sexratio_stats_MS_updated.xls – this is a spreadsheet with summarized statistics for sex ratio changes across the transgenerational experiment.
8.  Surv_RT_statistics_MS.xls - this is a spreadsheet with summarized statistics for survival in the reciprocal transplant part of the experiment.
9.	Survival_stats_MS_updated.xls – this is a spreadsheet with summarized statistics for survival to adulthood across the transgenerational experiment.


IV.	LFS
A.	This directory represents large files that exceed file upload limits for GitHub and are stored using Git LFS (Large File Storage).
1.	filtered_variants_Ahudsonica_MS.txt – this is a large, 300 Megabyte tab-delimited text file for evaluating Fst values. There are 48 columns. Details on column information can be found at https://varscan.sourceforge.net/using-varscan.html. To download the file, navigate to: https://github.com/dam-lab/hudsonica_transgenerational_MS/blob/master/LFS/filtered_variants_Ahudsonica_MS.txt. 
