<div align="center">
    <h1>Data Analysis Pipeline</h1>
    <h4>An R pipeline for technical and biological data assessment.</h4>
    <p>
        <img src="https://img.shields.io/badge/Submitted in:-Proteomics_--_Clinical_Applications-blue" alt="Submitted in Proteomics - Clinical Applications" />
    </p>
    <p>
        <img src="https://img.shields.io/badge/Submitted:-11.10.2014-blue" alt="Submission: 11.10.2024" />
        <img src="https://img.shields.io/badge/Revised:-03.05.2025-orange" alt="Revision: 03.05.2025" />
        <img src="https://img.shields.io/badge/Accepted:-08.07.2025-green" alt="Accepted: 08.07.2025" />
        <img src="https://img.shields.io/badge/Published:-pending-green" alt="Published: 08.07.2025" />
    </p>
</div>

<!-- Table of Contents -->
# Table of Contents

- [About the Project](#about-the-project)
  * [Overall Project Description](#overall-project-description)
  * [About the Present Experiment](#about-the-present-experiment)
- [Run the Code](#run-the-code)

<br><br>

<!-- About the Project -->
# About the Project

## Overall Project Description

Depression is a pressing global issue with a prevalence of 3.8&nbsp;% worldwide, leading to severe social and economic consequences. Due to the complexity of the disease, research into prevention, diagnosis, and treatment options faces significant challenges. Understanding the intricate pathways by which depression develops and progresses is critical to expanding our knowledge base and developing effective intervention strategies. One area of interest is the possible role of hormonal contraception as a risk factor. While some studies suggest a link between hormonal contraceptives and depression, others find no association, leading to controversy in the medical literature.

Our project takes a biochemical approach, focusing on the effects of steroid hormones on neuronal cells in vitro to address potential biases in previous research. By treating neural progenitor cells with various hormonal contraceptives and analyzing resultant protein changes using quantitative proteomics, we aim to uncover insights into the pharmacological effects of steroids on the neuronal proteome.

We hope to provide a significant contribution to clarifying the possible influence of steroids on the development of depression. This knowledge is of utmost importance due to the widespread and long-term use of oral contraceptives by healthy people and could lead to better risk-benefit assessments.

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

## About the Present Experiment

The following descriptions are only intended to give an overall impression of what this code is used for. For further information, please refer to the published article. *(Note: The article got accepted, but is not published yet. This file will be updated with a link to the article once it is published.)*

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

### Cell culture

- 2D cell culture
- Cell line: ReNcell&nbsp;VM
- Treatments: 100&nbsp;ng/mL of ethinyl estradiol, levonorgestrel, their combination, and S-23 in DMSO
- Controls: 100&nbsp;ppm DMSO & pure cell culture medium
- Treatment duration: 14 days

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

### Sample Preparation

- Aliquot 1 – label-based quantification (LBQ):
  - lysis
  - reduction
  - alkylation
  - digestion with Trypsin/LysC
  - labeling with isobaric mass tags: TMTpro 16plex
  - pooling
  - desalting and pre-fractionation at high pH

- Aliquot 2 – label-based quantification (LFQ):
  - lysis
  - reduction
  - alkylation
  - digestion with Trypsin

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

### nano LC-HR-MS/MS Analysis

- LBQ:
  - Institute of Biochemistry at the German Sport University Cologne, Cologne, Germany
  - Q&nbsp;Exactive
  - Top-15 DDA
  - Full scan: Res: 70,000&nbsp;FWHM (at _m_/_z_ 200); AGC: 3e6, IT: 50&nbsp;ms; _m_/_z_ 375–1500
  - MS/MS: Res: 35,000&nbsp;FWHM (at _m_/_z_ 200); AGC: 2e5; IT: 250&nbsp;ms; Isolation width: _m_/_z_ 0.7
  - 0&nbsp;% to 40&nbsp;% ACN in 180&nbsp;min

- LFQ:
  - Leibnitz-Institut für Analytische Wissenschaften—ISAS—e.V., Dortmund, Germany
  - Q&nbsp;Exactive&nbsp;HF
  - Top-15 DDA
  - Full scan: Res: 60,000&nbsp;FWHM (at _m_/_z_ 200); AGC: 3e6, IT: 120&nbsp;ms; _m_/_z_ 300–1500
  - MS/MS: Res: 15,000&nbsp;FWHM (at _m_/_z_ 200); AGC: 5e4; IT: 50&nbsp;ms; Isolation width: _m_/_z_ 1.6
  - 0&nbsp;% to 30&nbsp;% ACN in 120&nbsp;min

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

### Data Analysis

- Database search:
  - Proteome Discoverer 3.0 with Sequest&nbsp;HT
  - UniProtKB/Swiss-Prot human reference proteome (Tax. ID: 9606)
- Statistical analysis:
  - [`MSstats`](https://github.com/Vitek-Lab/MSstats) and [`MSststsTMT`](https://github.com/Vitek-Lab/MSstatsTMT)
  - [`clusterProfiler`](https://github.com/YuLab-SMU/clusterProfiler) and [`DOSE`](https://github.com/YuLab-SMU/DOSE)

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>

<!-- Running the code -->
# Run the code

To execute the pipeline yourself, follow these steps:

**Setup:** Ensure you have installed JupyterLab with an R kernel. For simplicity and to avoid compatibility issues, you can use [my Docker container](https://github.com/SamThilmany/JupyterLab-with-R_Docker-Environment).

**Prepare the Code:** Download or fork this repository. If you want, rename the folder (e.g., to `thilmany-etal`) and move it into the `notebooks/` directory of your JupyterLab Docker environment.

**Run the Pipeline:** Start the Docker container and access JupyterLab. Navigate to `notebooks/thilmany-etal/`, open the Jupyter Notebook (`Analysis.ipynb`), and select "Restart Kernel and Run All Cells..." from the "Run" menu.

Depending on Docker's resource allocation and your network connection, the process may take a while to finish, so do yourself a favor and have a cup of coffee. ☕️

<p align="right">(<a href="#table-of-contents">back to top</a>)</p>
