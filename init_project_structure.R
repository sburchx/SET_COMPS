# Load or install required package
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

# Define subfolders for each analysis type
modules <- c("flow", "rnaseq")
subdirs <- c("00_scripts", "01_data", "02_results", "03_rdata")

# Loop through and create folder structure
for (module in modules) {
  for (sub in subdirs) {
    dir_path <- here(module, sub)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      message("Created: ", dir_path)
    }
  }
}

# Create full README.md with template if it doesn't exist
readme_path <- here::here("README.md")

if (!file.exists(readme_path)) {
  readme_contents <- c(
    "# EXPERIMENT Project",
    "### Project details:",
    "- Cell Lines: ",
    "- Treatments:",
    "- Collection Days:",
    "",
    "This repository contains analysis for experiment EXPERIMENT, including:",
    "",
    "- Flow cytometry analysis (`flow/`)",
    "- Bulk RNA-seq analysis (`rnaseq/`)",
    "- ATAC-seq analysis (`atac/`)",
    "",
    "## 📁 Project Structure",
    "",
    "Each analysis folder contains:",
    "",
    "- `00_scripts/`: R scripts for preprocessing, analysis, and visualization",
    "- `01_data/`: Raw and processed input files",
    "- `02_results/`: Output plots, tables, and summary stats",
    "- `03_rdata/`: Saved R objects for reuse",
    "",
    "```",
    "EXPERIMENT/",
    "├── EXPERIMENT.Rproj",
    "├── README.md",
    "├── flow/",
    "│   ├── 00_scripts/",
    "│   ├── 01_data/",
    "│   ├── 02_results/",
    "│   ├── 03_rdata/",
    "├── rnaseq/",
    "│   ├── 00_scripts/",
    "│   ├── 01_data/",
    "│   ├── 02_results/",
    "│   ├── 03_rdata/",
    "├── atac/",
    "│   ├── 00_scripts/",
    "│   ├── 01_data/",
    "│   ├── 02_results/",
    "│   ├── 03_rdata/",
    "```",
    "",
    "#### 📦 Using the `here` Package for File Paths",
    "This project uses the [`here`](https://cran.r-project.org/package=here) package to manage file paths in a clean and reproducible way. Instead of using `setwd()` or hardcoding file paths, `here()` automatically sets the project root to the location of the `.Rproj` file.",
    "",
    "- Works consistently regardless of where your script is located",
    "- Avoids breaking paths when moving or sharing the project",
    "- Keeps code clean and portable",
    "",
    "#### 🔧 How to use it:",
    "",
    "1. **Install (only once):**",
    "```r",
    "install.packages(\"here\")",
    "```",
    "",
    "2. **Load it in any script or RMarkdown:**",
    "```r",
    "library(here)",
    "```",
    "",
    "3. **Use it for reading/writing files relative to the project root:**",
    "```r",
    "# Read data",
    "data <- read.csv(here(\"00_data\", \"counts.csv\"))",
    "",
    "# Save plot",
    "ggsave(here(\"02_results\", \"plot.png\"))",
    "```",
    "",
    "By using `here()`, you don’t need to use `setwd()` or worry about relative paths — everything will just work as long as you open the project using the `.Rproj` file.",
    "",
    "",
    "```",
    "",
    "## ✍️ Author",
    "**Shalise Burch** ",
    "PhD Candidate @ University of Chicago ",
    "Email: burchs@uchicago.com ",
    "GitHub: sburchx"
  )
  
  writeLines(readme_contents, readme_path)
  message("Created: README.md")
} else {
  message("README.md already exists")
}
