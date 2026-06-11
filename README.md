# T21 Project

### Experiment details:

-   [Cell Lines:]{.underline}
    -   **SET2**: WiCell - Bhattacharyya Collection

        -   Trisomy 21 hiPSCUWWC1-DS1

        -   Disomy 21 (Isogenic) hiPSC UWWC1-DS2U

        -   Established from a 1 year old male

        -   \*Isogenic hiPSC control in Set2 was achieved by spontaneous loss of trisomy.

    -   **SET3**: WiCell - Bhattacharyya Collection

        -   Trisomy 21 hiPSC WC-24-02-DS-M

        -   Disomy 21 (Isogenic) hiPSC WC-24-02-DS-B

        -   Established from a 25 year old female.

        -   \*Isogenic hiPSC control in Set3 was achieved by mosaic trisomy separation.

    -   **SET4**: Bhattacharyya Direct

        -   Trisomy 21 hiPSC WC68-35

        -   Disomy 21 (Isogenic) hiPSC WC68-18

        -   Established from a male.

        -   \*Isogenic hiPSC control in Set3 was achieved by mosaic trisomy separation.
-   [Treatments:]{.underline}
    -   No Treatment (CTRL)

    -   1uM Retinoic acid (RA) D3-D7
-   [Collection Days:]{.underline}
    -   Flow:
        -   SET2: SAB02, SAB06
            -   CTRL/RA: D3 -\> D10
        -   SET3: SAB01, SAB02, SAB04, SAB08
            -   CTRL: D3 -\> D10
            -   RA: D3 -\> D10
        -   SET4: SAB08
            -   CTRL: D3 -\> D10
            -   RA: D3 -\> D10
    -   RNAseq
        -   SET2:
            -   CTRL/RA: D3, D4, D5, D7
        -   SET3:
            -   CTRL/RA: D3, D4, D5, D6, D7, D8, D10
        -   SET4:
            -   CTRL/RA: D3, D4, D6, D8, D10

This repository contains analysis for multiple experiments, including:

-   Flow cytometry analysis (`flow/`)
-   Bulk RNA-seq analysis (`rnaseq/`)

## 📁 Project Structure

Each analysis folder contains:

-   `00_scripts/`: R scripts for preprocessing, analysis, and visualization
-   `01_data/`: Raw and processed input files
-   `02_results/`: Output plots, tables, and summary stats
-   `03_rdata/`: Saved R objects for reuse

```         
EXPERIMENT/
├── EXPERIMENT.Rproj
├── README.md
├── flow/
│   ├── 00_scripts/
│   ├── 01_data/
│   ├── 02_results/
│   ├── 03_rdata/
├── rnaseq/
│   ├── 00_scripts/
│   ├── 01_data/
│   ├── 02_results/
│   ├── 03_rdata/
```

This project uses the [`here`](https://cran.r-project.org/package=here) package to manage file paths in a clean and reproducible way. Instead of using `setwd()` or hardcoding file paths, `here()` automatically sets the project root to the location of the `.Rproj` file.

#### ✅ Why use `here()`?

-   Works consistently regardless of where your script is located
-   Avoids breaking paths when moving or sharing the project
-   Keeps code clean and portable

#### 🔧 How to use it:

1.  **Install (only once):**

``` r
install.packages("here")
```

2.  **Load it in any script or RMarkdown:**

``` r
library(here)
```

3.  **Use it for reading/writing files relative to the project root:**

``` r
# Read data
data <- read.csv(here("00_data", "counts.csv"))

# Save plot
ggsave(here("02_results", "plot.png"))
```

By using `here()`, you don’t need to use `setwd()` or worry about relative paths — everything will just work as long as you open the project using the `.Rproj` file.

## ✍️ Author

**Shalise Burch**\
PhD Candidate \@ University of Chicago

-   Email:[burchs\@uchicago.com](mailto:burchs@uchicago.com){.email}
-   GitHub: sburchx
