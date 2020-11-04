# Longitudinal Olink Proteomics Data and Code

Repository for the paper: “Longitudinal proteomic profiling of high-risk patients with COVID-19 reveals markers of severity and predictors of fatal disease” by Gisby et al, MedRxiv link to come. 

Here, we have stored:
* A public version of the dataset used in the paper (`data`)
* An RMD script replicating our major analyses (`scripts`)
* A modified version of the OlinkAnalyze package containing functions used to fit models and generate plots (`OlinkAnalyze_Modified`, https://github.com/Olink-Proteomics/OlinkRPackage/tree/master/OlinkAnalyze)

## Description of the Dataset

The dataset has been split into a series of four CSV files:
* `plasma_npx_level.csv` - Protein expression (NPX) data in long table format for the primary (plasma) cohort
* `plasma_sample_level.csv` - Sample level phenotypic data for the primary (plasma) cohort
* `serum_npx_level.csv` - Protein expression (NPX) data in long table format for the validation (serum) cohort
* `serum_sample_level.csv` - Sample level phenotypic data for the validation (serum) cohort

Note that one individual ("C105") of the primary cohort was sampled both before and after
their first postivie swab, so the first sample is marked as "NEGATIVE" for `Case_Control`.

### Protein Expression Data Features

For each of the primary and validation cohorts we provide the following columns:

Column Name | Data Type | Description
| :---: | :---: | :---:
SampleID  | Character | Unique identifier for samples
Individual_ID | Character | Unique identifier for individuals
UniProt | Character | UniProt ID for each protein
GeneID | Character | GeneID associated with each protein
Assay | Character | Name assigned to the protein by Olink
NPX | Numeric | Normalised protein expression (log2)
Panel | Character | The name of the panel each protein belongs to
Index | Integer | Unique sample identifier

### Sample Level Data Features

The features of the sample level datasets include:

Column Name | Data Type | Description
| :---: | :---: | :---:
SampleID  | Character | Unique identifier for samples
Individual_ID | Character | Unique identifier for individuals
Plate_ID | Character | The assay plate to which the sample was assigned (1-5)
Case_Control | Character | Whether the individual was COVID `POSITIVE` or `NEGATIVE` at time of sampling. We sampled one individual ("C105") prior to their first positive swab, so they have both `POSITIVE` and `NEGATIVE` time points.
WHO_Severity_Peak | Character | The peak (WHO) severity for the patient over the disease course
WHO_Severity_Contemporaneous | Character | The (WHO) severity at time of sampling
Sex | Character | The individual's sex (M or F)
Ethnicity | Character | The individual's ethnicity (White, Black, South Asian, Asian (other) or Other)
Age | Character | Age in 20 year bins
Ever_Admitted | Logical | Whether the individual was admitted to hospital due to COVID at any point
Fatal_Disease | Logical | Whether the disease was fatal
Time_From_First_Symptoms | Integer | The number of days since the individual first experienced COVID symptoms at time of sampling
Time_Of_Death_From_First_Symptoms | Integer | The number of days between the individual's first experienced COVID symptoms and time of death
Time_From_First_Swab | Integer | The number of days since the individual's first positive swab was taken at time of sampling
Time_Of_Death_From_First_Swab | Integer | The number of days between the individual's first positive swab and time of death
WCC | Numeric | Individual's blood white cell count (x10^9/L) at time of sampling
Neutrophil | Numeric | Individual's blood neutrophil count (x10^9/L) at time of sampling
Monocyte | Numeric | Individual's blood monocyte count (x10^9/L) at time of sampling
Lymophocyte | Numeric | Individual's blood lymphocyte count (x10^9/L) at time of sampling
CRP | Numeric | Individual's blood C-Reactive Protein concentration (mg/L) at time of sampling
DDimer | Integer | Individual's blood D-Dimer concentration (ug/L) at time of sampling
Ferritin | Integer | Individual's blood ferritin concentration (ug/L) at time of sampling
Troponin | Integer | Individual's blood troponin concentration (ng/L) at time of sampling

## Replicating the Analyses

We have provided data and R code for replicating the major analyses. The RMD 
script can be run (`scripts/example_script.Rmd`) to re-create the main results.
Note that we used exact age, rather than binned age, in our models so results may 
differ slightly.

The RMD script uses a modified OlinkAnalyze script to run the analyses and this
must be installed prior to running. Code for this is also available in this 
repository if of interest. 
