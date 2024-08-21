---
output:
  html_document:
    theme: united
---
#### **NanoTAXI App allows users to process Nanopore 16S Sequencing Data in real-time and visualize**

- **Explore** the app's features with the example count data set pre-loaded by clicking on the Input Data tab.
- **Upload** your own samples list in the "Input Data" tab.


## <a name="Features"></a> **Features**

<!-- <style>
  .features {
  display: grid;
  grid-template-columns: repeat(2, 1fr);
  row-gap: 10px;
  column-gap: 10px;
  }
</style> -->


<img src="Drug_class.png" alt="Drug Class vs. ARG Abundance" style="width:45%"/>
<img src="Resistance_mechanism.png" alt="Resistance Mechanism" style="width:45%"/>
<img src="ARG_Cohort_Richness.png" alt="Cohort ARG Richness" style="width:45%"/>
<img src="ARG_Cohort_Abundance.png" alt="Cohort ARG Abundance" style="width:45%"/>
<img src="Alpha_diversity.png" alt="Alpha Diversity" style="width:45%"/>
<img src="Abundance_Kruskal_Walis.png" alt="Abundance Comparison" style="width:45%"/>
<img src="HeatMap.png" alt="HeatMap" style="width:45%"/>
<img src="PCA_plot.png" alt="PCA Plot" style="width:45%"/>


## **Data Visualization**

- ARG Mechanisms (Drug Class vs ARG Abundance, Resistance Mechanisms, ARG Richness per Bacterial Species, ARG Abundance per Bacterial Species)
- Alpha Diversity and Abundance Comparison (Box-plots)
- Clustering (HeatMap, PCA Plot)

## <a name="DataFormat"></a> **Data Format**

- Must be a .csv **comma-separated-value** file.
- File must have two columns with the headers 1: Sample_Id, 2: Group.
- 1: Sample_Id: The name of the sample.
- 2: Group: The group that the sample belongs to.
- The first row of the file is the header.

## <a name="InputData"></a> **Input Data**
- Each row from the first column representes the sample name.
- Each row from the second column representes the group (Control or Case).

<style>
  .sample_info {
    height: 200px;
    overflow-y: scroll;
    border: 1px solid #000000; 
    width: 55%;
    margin-left:5%;
    margin-right:5%
  }

  .sample_info th, .sample_info td {
    border: 1px solid #000000;
    border-spacing: 0;
    text-align: center;
    padding-top: 10px; padding-left:70px; padding-right: 70px; padding-bottom: 10px
  }
</style>

<div class="sample_info">

| Sample_Id            | Group    
|-----------------|--------
| SRR5298196        |   Case   
| SRR5298197        |   Case   
| SRR5298198        |   Case   
| SRR5298199        |   Case   
| SRR5298200        |   Case   
| SRR5298201        |   Case   
| SRR5298202        |   Case   
| SRR5298203        |   Case   
| SRR5298204        |   Case   
| SRR5298205        |   Case   
| SRR5298206        |   Case   
| SRR5298207        |   Case   
| SRR5298208        |   Case   
| SRR5298209        |   Case   
| SRR5298210        |   Case   
| SRR5298211        |   Case   
| SRR5298212        |   Case   
| SRR5298213        |   Case   
| SRR5298214        |   Control
| SRR5298215        |   Control
| SRR5298216        |   Control
| SRR5298217        |   Control
| SRR5298218        |   Control
| SRR5298219        |   Control
| SRR5298220        |   Control
| SRR5298221        |   Control
| SRR5298222        |   Control
</div>


## <a name="OutputData"></a> **Output Data**
- Each row of the First column represents the ARG term.
- Additional columns provide information about 
1) Drug Class 
2) Resistance Mechanism 
3) AMR Gene Family
4) Bacterial Classification, etc.
5) Normalized Counts


<style>
  .output_data {
    height: 400px;
    overflow-y: scroll;
    border: 1px solid #cccccc;
  }
  .output_data th, .output_data td {
      border: 1px solid #000000;
      border-spacing: 0;
      text-align: center;
      padding: 8px;
      }
</style>
<div class="output_data">

|**ARO_term**|**ARG_length**|**Counts**|**Percentage_Identity**|**Drug_Class**|**Resistance_Mechanism**| **AMR_Gene_Family**|**Percentage_Coverage**|**Classification**|**Sample_Id**|**Normalized_counts**|**Family**|**Group**|
|:--------------------------------------------------------------------------------:|:----------------:|:------------:|:-------------------------:|:--------------------------------------------:|:------------------------------------:|:---------------------------------------------------------------------------------------:|:-------------------------:|:----------------------------------:|:---------------:|:-----------------------:|:---------------------------:|:-----------:|
| "ErmX"                                                                           | 855              | 160          | 90.14                     | "macrolide antibiotic"                       | "antibiotic target alteration"       | "Erm 23S ribosomal RNA methyltransferase"                                               | 100                       | "Rothia"                           | "SRR5298196"    | 1333.29                 | "Micrococcaceae"            | "Case"      |
| "AAC(3)-IIc"                                                                     | 801              | 23           | 96.24                     | "aminoglycoside antibiotic"                  | "antibiotic inactivation"            | "AAC(3)"                                                                                | 93.01                     | "Rothia"                           | "SRR5298196"    | 204.58                  | "Micrococcaceae"            | "Case"      |
| "tetM"                                                                           | 1920             | 8268         | 98.9                      | "tetracycline antibiotic"                    | "antibiotic target protection"       | "tetracycline-resistant ribosomal protection protein"                                   | 100                       | "Streptococcus sp000187445"        | "SRR5298196"    | 30680.97                | "Streptococcaceae"          | "Case"      |
| "msrE"                                                                           | 1476             | 5539         | 100                       | "macrolide antibiotic"                       | "antibiotic target protection"       | "ABC-F ATP-binding cassette ribosomal protection protein"                               | 100                       | "Streptococcus sp000187445"        | "SRR5298196"    | 26737.13                | "Streptococcaceae"          | "Case"      |
| "mphE"                                                                           | 885              | 3190         | 100                       | "macrolide antibiotic"                       | "antibiotic inactivation"            | "macrolide phosphotransferase (MPH)"                                                    | 100                       | "Streptococcus sp000187445"        | "SRR5298196"    | 25681.31                | "Streptococcaceae"          | "Case"      |
| "ugd"                                                                            | 1167             | 3127         | 99.23                     | "peptide antibiotic"                         | "antibiotic target alteration"       | "pmr phosphoethanolamine transferase"                                                   | 100                       | "unclassified"                     | "SRR5298196"    | 19090.92                | "unclassified"              | "Case"      |
| "gadW"                                                                           | 729              | 3585         | 99.59                     | "macrolide antibiotic"                       | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 35037.36                | "unclassified"              | "Case"      |
| "mdtA"                                                                           | 1248             | 2529         | 99.04                     | "aminocoumarin antibiotic"                   | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 14437.9                 | "unclassified"              | "Case"      |
| "mdtB"                                                                           | 3123             | 5630         | 99.81                     | "aminocoumarin antibiotic"                   | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 12844.17                | "unclassified"              | "Case"      |
| "mdtC"                                                                           | 3060             | 5628         | 99.61                     | "aminocoumarin antibiotic"                   | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 99.41                     | "unclassified"                     | "SRR5298196"    | 13103.96                | "unclassified"              | "Case"      |
| "baeS"                                                                           | 1404             | 2396         | 99.36                     | "aminoglycoside antibiotic"                  | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 12158.77                | "unclassified"              | "Case"      |
| "baeR"                                                                           | 723              | 1523         | 98.75                     | "aminoglycoside antibiotic"                  | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 15008.3                 | "unclassified"              | "Case"      |
| "acrB"                                                                           | 3150             | 8986         | 100                       | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 20324.77                | "unclassified"              | "Case"      |
| "Escherichia coli acrA"                                                          | 1194             | 4109         | 99.75                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 24518.93                | "unclassified"              | "Case"      |
| "evgS"                                                                           | 3594             | 2097         | 97.41                     | "macrolide antibiotic"                       | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 4157.1                  | "unclassified"              | "Case"      |
| "mdtG"                                                                           | 1227             | 3028         | 99.75                     | "fosfomycin"                                 | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 17582.52                | "unclassified"              | "Case"      |
| "mdtH"                                                                           | 1209             | 2645         | 99.75                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 15587.24                | "unclassified"              | "Case"      |
| "PmrF"                                                                           | 969              | 1127         | 99.69                     | "peptide antibiotic"                         | "antibiotic target alteration"       | "pmr phosphoethanolamine transferase"                                                   | 100                       | "unclassified"                     | "SRR5298196"    | 8286.48                 | "unclassified"              | "Case"      |
| "Escherichia coli ampH beta-lactamase"                                           | 1158             | 3572         | 99.74                     | "cephalosporin"                              | "antibiotic inactivation"            | "ampC-type beta-lactamase"                                                              | 100                       | "unclassified"                     | "SRR5298196"    | 21977.22                | "unclassified"              | "Case"      |
| "emrB"                                                                           | 1539             | 522          | 99.8                      | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 2416.58                 | "unclassified"              | "Case"      |
| "Escherichia coli ampC1 beta-lactamase"                                          | 1305             | 1870         | 100                       | "cephalosporin"                              | "antibiotic inactivation"            | "ampC-type beta-lactamase"                                                              | 100                       | "unclassified"                     | "SRR5298196"    | 10209.42                | "unclassified"              | "Case"      |
| "kdpE"                                                                           | 678              | 892          | 100                       | "aminoglycoside antibiotic"                  | "antibiotic efflux"                  | "kdpDE"                                                                                 | 100                       | "unclassified"                     | "SRR5298196"    | 9373.57                 | "unclassified"              | "Case"      |
| "mdtN"                                                                           | 1032             | 3330         | 98.83                     | "nucleoside antibiotic"                      | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 22989.75                | "unclassified"              | "Case"      |
| "CRP"                                                                            | 633              | 2845         | 99.52                     | "macrolide antibiotic"                       | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 32021.99                | "unclassified"              | "Case"      |
| "YojI"                                                                           | 1644             | 3309         | 99.63                     | "peptide antibiotic"                         | "antibiotic efflux"                  | "ATP-binding cassette (ABC) antibiotic efflux pump"                                     | 100                       | "unclassified"                     | "SRR5298196"    | 14340.51                | "unclassified"              | "Case"      |
| "eptA"                                                                           | 1644             | 4729         | 100                       | "peptide antibiotic"                         | "antibiotic target alteration"       | "pmr phosphoethanolamine transferase"                                                   | 100                       | "unclassified"                     | "SRR5298196"    | 20494.5                 | "unclassified"              | "Case"      |
| "marA"                                                                           | 384              | 907          | 99.21                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 16828.52                | "unclassified"              | "Case"      |
| "bacA"                                                                           | 822              | 2541         | 99.63                     | "peptide antibiotic"                         | "antibiotic target alteration"       | "undecaprenyl pyrophosphate related proteins"                                           | 100                       | "unclassified"                     | "SRR5298196"    | 22024.33                | "unclassified"              | "Case"      |
| "gadX"                                                                           | 825              | 3287         | 97.45                     | "macrolide antibiotic"                       | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 28386.74                | "unclassified"              | "Case"      |
| "eptB"                                                                           | 1692             | 1472         | 90.07                     | "peptide antibiotic"                         | "antibiotic target alteration"       | "pmr phosphoethanolamine transferase"                                                   | 98.08                     | "unclassified"                     | "SRR5298196"    | 6198.37                 | "unclassified"              | "Case"      |
| "QnrS1"                                                                          | 657              | 1648         | 100                       | "fluoroquinolone antibiotic"                 | "antibiotic target protection"       | "quinolone resistance protein (qnr)"                                                    | 100                       | "unclassified"                     | "SRR5298196"    | 17871.53                | "unclassified"              | "Case"      |
| "ArnT"                                                                           | 1656             | 277          | 89.29                     | "peptide antibiotic"                         | "antibiotic target alteration"       | "pmr phosphoethanolamine transferase"                                                   | 100                       | "unclassified"                     | "SRR5298196"    | 1191.76                 | "unclassified"              | "Case"      |
| "AcrF"                                                                           | 3105             | 11556        | 99.52                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 26516.47                | "unclassified"              | "Case"      |
| "QnrB4"                                                                          | 648              | 1922         | 100                       | "fluoroquinolone antibiotic"                 | "antibiotic target protection"       | "quinolone resistance protein (qnr)"                                                    | 100                       | "unclassified"                     | "SRR5298196"    | 21132.37                | "unclassified"              | "Case"      |
| "DHA-1"                                                                          | 1140             | 3186         | 100                       | "cephalosporin"                              | "antibiotic inactivation"            | "DHA beta-lactamase"                                                                    | 100                       | "unclassified"                     | "SRR5298196"    | 19911.81                | "unclassified"              | "Case"      |
| "sul1"                                                                           | 840              | 1151         | 100                       | "sulfonamide antibiotic"                     | "antibiotic target replacement"      | "sulfonamide resistant sul"                                                             | 100                       | "unclassified"                     | "SRR5298196"    | 9762.61                 | "unclassified"              | "Case"      |
| "sul2"                                                                           | 816              | 1230         | 100                       | "sulfonamide antibiotic"                     | "antibiotic target replacement"      | "sulfonamide resistant sul"                                                             | 100                       | "unclassified"                     | "SRR5298196"    | 10739.52                | "unclassified"              | "Case"      |
| "APH(3'')-Ib"                                                                    | 804              | 1609         | 99.63                     | "aminoglycoside antibiotic"                  | "antibiotic inactivation"            | "APH(3'')"                                                                              | 100                       | "unclassified"                     | "SRR5298196"    | 14258.37                | "unclassified"              | "Case"      |
| "APH(6)-Id"                                                                      | 837              | 2039         | 99.64                     | "aminoglycoside antibiotic"                  | "antibiotic inactivation"            | "APH(6)"                                                                                | 100                       | "unclassified"                     | "SRR5298196"    | 17356.48                | "unclassified"              | "Case"      |
| "tet(59)"                                                                        | 1203             | 476          | 99.25                     | "tetracycline antibiotic"                    | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 2819.1                  | "unclassified"              | "Case"      |
| "TolC"                                                                           | 1482             | 5186         | 99.8                      | "macrolide antibiotic"                       | "antibiotic efflux"                  | "ATP-binding cassette (ABC) antibiotic efflux pump"                                     | 99.6                      | "unclassified"                     | "SRR5298196"    | 24931.83                | "unclassified"              | "Case"      |
| "rsmA"                                                                           | 186              | 847          | 85.25                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 32444.44                | "unclassified"              | "Case"      |
| "SHV-66"                                                                         | 861              | 92           | 99.65                     | "carbapenem"                                 | "antibiotic inactivation"            | "SHV beta-lactamase"                                                                    | 100                       | "unclassified"                     | "SRR5298196"    | 761.3                   | "unclassified"              | "Case"      |
| "msbA"                                                                           | 1749             | 126          | 91.92                     | "nitroimidazole antibiotic"                  | "antibiotic efflux"                  | "ATP-binding cassette (ABC) antibiotic efflux pump"                                     | 100                       | "unclassified"                     | "SRR5298196"    | 513.28                  | "unclassified"              | "Case"      |
| "CTX-M-15"                                                                       | 876              | 78           | 100                       | "cephalosporin"                              | "antibiotic inactivation"            | "CTX-M beta-lactamase"                                                                  | 100                       | "unclassified"                     | "SRR5298196"    | 634.4                   | "unclassified"              | "Case"      |
| "Klebsiella pneumoniae KpnH"                                                     | 1539             | 136          | 93.23                     | "macrolide antibiotic"                       | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 629.61                  | "unclassified"              | "Case"      |
| "dfrA14"                                                                         | 459              | 971          | 100                       | "diaminopyrimidine antibiotic"               | "antibiotic target replacement"      | "trimethoprim resistant dihydrofolate reductase dfr"                                    | 96.82                     | "unclassified"                     | "SRR5298196"    | 15072.19                | "unclassified"              | "Case"      |
| "tet(A)"                                                                         | 1275             | 43           | 98.57                     | "tetracycline antibiotic"                    | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 240.29                  | "unclassified"              | "Case"      |
| "H-NS"                                                                           | 414              | 1260         | 100                       | "macrolide antibiotic"                       | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 21684.03                | "unclassified"              | "Case"      |
| "TEM-1"                                                                          | 861              | 4885         | 100                       | "monobactam"                                 | "antibiotic inactivation"            | "TEM beta-lactamase"                                                                    | 100                       | "unclassified"                     | "SRR5298196"    | 40423.25                | "unclassified"              | "Case"      |
| "emrR"                                                                           | 531              | 53           | 92.57                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "major facilitator superfamily (MFS) antibiotic efflux pump"                            | 100                       | "unclassified"                     | "SRR5298196"    | 711.13                  | "unclassified"              | "Case"      |
| "AAC(6')-Ig"                                                                     | 438              | 20           | 100                       | "aminoglycoside antibiotic"                  | "antibiotic inactivation"            | "AAC(6')"                                                                               | 100                       | "unclassified"                     | "SRR5298196"    | 325.33                  | "unclassified"              | "Case"      |
| "Enterobacter cloacae acrA"                                                      | 1194             | 121          | 96.22                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 722.02                  | "unclassified"              | "Case"      |
| "ErmB"                                                                           | 738              | 2431         | 97.96                     | "macrolide antibiotic"                       | "antibiotic target alteration"       | "Erm 23S ribosomal RNA methyltransferase"                                               | 98.79                     | "unclassified"                     | "SRR5298196"    | 23469.2                 | "unclassified"              | "Case"      |
| "ramA"                                                                           | 357              | 27           | 99.16                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 95.97                     | "unclassified"                     | "SRR5298196"    | 538.85                  | "unclassified"              | "Case"      |
| "Escherichia coli EF-Tu mutants conferring resistance to Pulvomycin"             | 1185             | 10119        | 99.75                     | "elfamycin antibiotic"                       | "antibiotic target alteration"       | "elfamycin resistant EF-Tu"                                                             | 96.33                     | "unclassified"                     | "SRR5298196"    | 60839.97                | "unclassified"              | "Case"      |
| "Escherichia coli UhpT with mutation conferring resistance to fosfomycin"        | 1392             | 4350         | 95.46                     | "fosfomycin"                                 | "antibiotic target alteration"       | "antibiotic-resistant UhpT"                                                             | 100                       | "unclassified"                     | "SRR5298196"    | 22264.85                | "unclassified"              | "Case"      |
| "Escherichia coli acrR with mutation conferring multidrug antibiotic resistance" | 648              | 2200         | 100                       | "fluoroquinolone antibiotic"                 | "antibiotic target alteration"       | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 24188.97                | "unclassified"              | "Case"      |
| "Escherichia coli marR mutant conferring antibiotic resistance"                  | 435              | 1112         | 98.61                     | "fluoroquinolone antibiotic"                 | "antibiotic target alteration"       | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "unclassified"                     | "SRR5298196"    | 18213.16                | "unclassified"              | "Case"      |
| "Escherichia coli soxR with mutation conferring antibiotic resistance"           | 465              | 1719         | 98.7                      | "fluoroquinolone antibiotic"                 | "antibiotic target alteration"       | "ATP-binding cassette (ABC) antibiotic efflux pump"                                     | 100                       | "unclassified"                     | "SRR5298196"    | 26338.6                 | "unclassified"              | "Case"      |
| "Escherichia coli soxS with mutation conferring antibiotic resistance"           | 324              | 1083         | 100                       | "fluoroquinolone antibiotic"                 | "antibiotic target alteration"       | "ATP-binding cassette (ABC) antibiotic efflux pump"                                     | 100                       | "unclassified"                     | "SRR5298196"    | 23815.14                | "unclassified"              | "Case"      |
| "CRP"                                                                            | 633              | 363          | 95.24                     | "macrolide antibiotic"                       | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "Vibrio cholerae"                  | "SRR5298196"    | 4085.76                 | "Vibrionaceae"              | "Case"      |
| "Vibrio cholerae varG"                                                           | 1125             | 365          | 100                       | "carbapenem"                                 | "antibiotic inactivation"            | "subclass B1 Vibrio cholerae varG beta-lactamase"                                       | 95.9                      | "Vibrio cholerae"                  | "SRR5298196"    | 2311.59                 | "Vibrionaceae"              | "Case"      |
| "almG"                                                                           | 822              | 233          | 100                       | "peptide antibiotic"                         | "antibiotic target alteration"       | "lipid A acyltransferase"                                                               | 100                       | "Vibrio cholerae"                  | "SRR5298196"    | 2019.55                 | "Vibrionaceae"              | "Case"      |
| "dfrA1"                                                                          | 474              | 288          | 99.36                     | "diaminopyrimidine antibiotic"               | "antibiotic target replacement"      | "trimethoprim resistant dihydrofolate reductase dfr"                                    | 100                       | "Vibrio cholerae"                  | "SRR5298196"    | 4328.96                 | "Vibrionaceae"              | "Case"      |
| "catB9"                                                                          | 630              | 256          | 100                       | "phenicol antibiotic"                        | "antibiotic inactivation"            | "chloramphenicol acetyltransferase (CAT)"                                               | 100                       | "Vibrio cholerae"                  | "SRR5298196"    | 2895.14                 | "Vibrionaceae"              | "Case"      |
| "tetO"                                                                           | 1920             | 4564         | 99.22                     | "tetracycline antibiotic"                    | "antibiotic target protection"       | "tetracycline-resistant ribosomal protection protein"                                   | 100                       | "Bifidobacterium longum"           | "SRR5298197"    | 27107.66                | "Bifidobacteriaceae"        | "Case"      |
| "dfrA14"                                                                         | 474              | 1018         | 100                       | "diaminopyrimidine antibiotic"               | "antibiotic target replacement"      | "trimethoprim resistant dihydrofolate reductase dfr"                                    | 100                       | "Bifidobacterium longum"           | "SRR5298197"    | 24491.59                | "Bifidobacteriaceae"        | "Case"      |
| "sul2"                                                                           | 816              | 976          | 100                       | "sulfonamide antibiotic"                     | "antibiotic target replacement"      | "sulfonamide resistant sul"                                                             | 100                       | "Bifidobacterium longum"           | "SRR5298197"    | 13639.78                | "Bifidobacteriaceae"        | "Case"      |
| "APH(6)-Id"                                                                      | 837              | 1606         | 99.28                     | "aminoglycoside antibiotic"                  | "antibiotic inactivation"            | "APH(6)"                                                                                | 100                       | "Bifidobacterium longum"           | "SRR5298197"    | 21881.03                | "Bifidobacteriaceae"        | "Case"      |
| "TEM-1"                                                                          | 861              | 1800         | 100                       | "monobactam"                                 | "antibiotic inactivation"            | "TEM beta-lactamase"                                                                    | 100                       | "Bifidobacterium longum"           | "SRR5298197"    | 23840.59                | "Bifidobacteriaceae"        | "Case"      |
| "marA"                                                                           | 384              | 399          | 100                       | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "Escherichia flexneri"             | "SRR5298197"    | 11849.21                | "Enterobacteriaceae"        | "Case"      |
| "rsmA"                                                                           | 186              | 343          | 85.25                     | "fluoroquinolone antibiotic"                 | "antibiotic efflux"                  | "resistance-nodulation-cell division (RND) antibiotic efflux pump"                      | 100                       | "Escherichia flexneri"             | "SRR5298197"    | 21029.49                | "Enterobacteriaceae"        | "Case"      |
</div>

## <a name="help"></a> **Additional Information**

Additional Information about the pipeline and help with all the visualization plots are provided under the "Additional Information" tab.

## **App Info**

The NanoTAXI App has been develoed by Nirmal Singh Mahar and Ishaan Gupta*.


#### **Please cite our App:**

[]

The source code of MetaShiny is available on [Github](https://github.com/Nirmal2310/NanoTAXI).

We would appreciate reports of any issues with the app via [Github](https://github.com/Nirmal2310/NanoTAXI/issues)