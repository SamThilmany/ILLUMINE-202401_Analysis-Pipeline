<div align="center">
    <h1>Data Analysis Pipeline</h1>
    <h4>An R pipeline for technical and biological data assessment.</h4>
    <p>
        <img src="https://img.shields.io/badge/Submitted in:-Proteomics_--_Clinical_Applications-blue" alt="Submitted in Proteomics - Clinical Applications" />
    </p>
    <p>
        <img src="https://img.shields.io/badge/Submitted:-11.10.2014-blue" alt="Submission: 11.10.2024" />
        <img src="https://img.shields.io/badge/Revised:-pending-orange" alt="Revision pending" />
        <img src="https://img.shields.io/badge/Accepted:-pending-green" alt="Acceptance pending" />
    </p>
</div>

<!-- Table of Contents -->
# Table of Contents

- [About the Project](#about-the-project)
  * [Overall Project Description](#overall-project-description)
  * [About the Present Experiment](#about-the-present-experiment)
- [Run the Code](#run-the-code)
  * [General Remarks on the Repository](#general-remarks-on-the-repository)
  * [Execute the Code Yourself](#execute-the-code-yourself)

<br><br>

<!-- About the Project -->
# About the Project

## Overall Project Description

Depression is a pressing global issue with a prevalence of 3.8&nbsp;% worldwide, leading to severe social and economic consequences. Due to the complexity of the disease, research into prevention, diagnosis, and treatment options faces significant challenges. Understanding the intricate pathways by which depression develops and progresses is critical to expanding our knowledge base and developing effective intervention strategies. One area of interest is the possible role of hormonal contraception as a risk factor. While some studies suggest a link between hormonal contraceptives and depression, others find no association, leading to controversy in the medical literature.

Our project takes a biochemical approach, focusing on the effects of steroid hormones on neuronal cells in vitro to address potential biases in previous research. By treating neural progenitor cells with various hormonal contraceptives and analyzing resultant protein changes using quantitative proteomics, we aim to uncover insights into the pharmacological effects of steroids on the neuronal proteome.

We hope to provide a significant contribution to clarifying the possible influence of steroids on the development of depression. This knowledge is of utmost importance due to the widespread and long-term use of oral contraceptives by healthy people and could lead to better risk-benefit assessments.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

## About the Present Experiment

The following descriptions are only intended to give an overall impression of what this code is used for. For further information, please refer to the published article. *(Note: The article will be submitted shortly. This file will be updated with a link to the article once it is published.)*

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

### Cell culture

In this experiment, the neural progenitor cell line ReNcell VM was cultured in a 2D cell culture and treated with various substances: ethynyl estradiol, levonorgestrel, the combination of ethynyl estradiol and levonorgestrel commonly found in oral contraceptive pills, and S-23, a drug candidate for the male contraceptive pill. The drugs were dissolved in DMSO and added to the cell culture medium at a final concentration of 100&nbsp;ng/mL. The final DMSO concentration in the cell culture medium was 100&nbsp;ppm.

One batch of cells was incubated with cell culture medium with added 100&nbsp;ppm DMSO and another batch was incubated with only cell culture medium. These cell cultures served as controls in the data analysis.

The cell culture medium was exchanged every other day, and the cells were split at approx 90&nbsp;% confluency. After 14 days, the cells were harvested, split into aliquots of $1 \times 10^6$ cells, and stored at -80&nbsp;°C.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

### Sample Preparation

One aliquot of each sample was processed for label-based quantification (LBQ) using isobaric mass tags (TMTpro 16plex), whereas another aliquot of each sample was processed for label-free quantification (LFQ).

In either case, the proteins were digested, reduced, and alkylated for bottom-up proteomics.

In the case of label-based quantification, the pooled sample was pre-fractionated at high pH into 8 fractions.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

### nano-LC-MS Analysis

For LC-MS data analysis, the samples were separated using nano-flow liquid chromatography and injected into a quadrupole-Orbitrap mass spectrometer operated in positive mode with a Top-15 data-dependent acquisition (DDA) method.

The LBQ samples were analyzed using a Thermo Scientific Q Exactive mass spectrometer at the Institute of Biochemistry at the German Sport University Cologne, Cologne, Germany. The LFQ samples were analyzed using a Thermo Scientific Q Exactive HF mass spectrometer at the Leibnitz-Institut für Analytische Wissenschaften—ISAS—e.V., Dortmund, Germany.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

### Data Analysis

The mass spectrometer's raw data was processed with Thermo Scientific Proteome Discoverer 3.0 (LBQ) and 3.1 (LFQ), using the Sequest HT and CHIMERYS search engines. The "Proteins", "Peptide Groups", "PSMs" and "MS/MS Spectrum Info" tables were exported as text files and further analyzed and visualized by the analysis pipeline in this repository. The data analysis involved multiple stages:
1. A technical evaluation was conducted for the LFQ and LBQ data to assess their quality.
2. The regulated proteins were biologically evaluated, including the search for gene-to-disease associations and enrichment analyses for GO (Gene Ontology) terms and KEGG (Kyoto Encyclopedia of Genes and Genomes) pathways.
3. The regulated proteins from all data sets were combined to obtain a more comprehensive view, and this combined data set was also subjected to biological analysis.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

<br><br>

<!-- Running the code -->
# Run the code

## General Remarks on the Repository

The code is organized with object-oriented programming principles in mind. In the `classes/` directory, you will find separate R files with S4 classes for a specific assessment, *e.g.*, the assessment of the AGC's fill percentage (`agc_fill.R`). To keep the code of these classes readable, the classes do not contain any code that manipulates data but only call functions responsible for those tasks. These functions can be found in the `utils/` directory with the `_helpers` suffix, *e.g.*, `agc_fill_helpers.R`.

Finally, the `notebooks/` directory contains three Jupyter Notebooks that combine all the necessary data analysis and visualization steps. The resulting files and images are saved in the `results/` directory.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

## Execute the Code Yourself

To execute the pipeline yourself, follow these steps:

**Setup:** Ensure you have installed JupyterLab with an R kernel. For simplicity and to avoid compatibility issues, you can use [my Docker container](https://github.com/SamThilmany/JupyterLab-with-R--Docker-Environment).

**Prepare the Code:** Download the code from this repository. Rename the folder (e.g., to `thilmany-etal`) and move it into the `notebooks/` directory of your JupyterLab Docker environment.

**Download Data:** Obtain the raw data from ProteomeXchange and place it into the appropriate folders within the `data/` directory. Ensure the file names match those specified in `utils/data.R` to avoid data loading errors.

**DISGENET API:** The gene-to-disease association analysis uses the [DISGENET](https://www.disgenet.com) gene-disease association network. Subscribe to a plan that provides access to the REST API and the R package (we used the Academic license). Copy your API key into the `.Renviron.example` file and rename it to `.Renviron`.

**Run the Pipeline:** Start the Docker container and access JupyterLab. Navigate to `notebooks/thilmany-etal/notebooks/`, open a Jupyter Notebook (e.g., `LBQ-analysis.ipynb`), and select "Restart Kernel and Run All Cells..." from the "Run" menu.

Depending on Docker's resource allocation, the process may take several minutes.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>
