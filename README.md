# Swift-BAT

Swift-BAT is a repository containing the analysis pipeline for searching and placing fluence limits on associated $\gamma$-ray emission from the Canadian Hydrogen Intensity Mapping Experiment/Fast Radio Burst (CHIME/FRB) telescope. 


If the corresponding SWIFT/BAT data does not exist locally, the pipeline will automatically download it.

The pipeline uses a boxcar search convolution over user-set time scales and energy bands around a user-set search window (centered on the CHIME trigger time). To search for FRBs. 

To produce fluence limits. The analysis pipeline produces a noise image, and a background subtracted count image for the search window. It then using the KBN method in astropy to determine a countlimit at 95%. It generates a SWIFT/BAT RSP file, and fits using XSPEC.

## Outputs

The pipeline produces two main outputs:

### Diagnostic plot

A diagnostic plot containing 5-6 figures per target, including:

- A signal-to-noise ratio (SNR) lightcurve
- A zoomed-in SNR lightcurve for the search window
- A peak SNR vs time scale plot
- A sky image (if one is produced)
- A histogram of the photon energy for the data
- A histogram of the photon energy for the search window

### Outcatalog JSON file

An `outcatalog.json` file containing the results of the analysis.

## Requirements

Swift-BAT requires the following Python packages:

- NumPy
- Matplotlib
- Astropy
- Re
- swifttools
- subprocess
- tqdm

In addition, Swift-BAT requires the following software packages:

- HEASoft
- Xspec

## Getting Started

To use Swift-BAT, you can clone the repository and run the pipeline on your local machine. Before running the pipeline, make sure you have installed the required Python and software packages listed above.

To run the pipeline, simply execute the `swift-bat.py` script with the appropriate command-line arguments.

A basic usage is

## Acknowledgments

This pipeline was developed by **Maxwell A. Fine** as part of AST425 (Undergraduate Thesis Class) at the University of Toronto. I was co-supervised by Ziggy and Paul, as well as Professor Bryan Gaensler. Aaron Tohuvavohu was instrumental in understanding HEASoft and working with BAT data. We acknowledge the use of the CHIME/FRB and SWIFT/BAT data in this analysis.

# Contact

If you have any questions or feedback, please contact **Maxwell Fine**:

- Website: https://afinemax.github.io/afinemax1/
- GitHub: https://github.com/afinemax
- Email: maxwell.fine@mail.utoronto.ca
