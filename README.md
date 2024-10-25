# Multivariate Signal Feature Extraction and Analysis

This project focuses on extracting and analyzing multivariate signal features from segmented data, including time-series statistics and signal characteristics. The extracted features are then visualized and compared across different data sets, enabling efficient classification and statistical analysis.

## Table of Contents

1. [Project Overview](#project-overview)
2. [Features Extracted](#features-extracted)
3. [Code Structure](#code-structure)
4. [Requirements](#requirements)
5. [Usage](#usage)
6. [Results and Visualizations](#results-and-visualizations)

## Project Overview

This project automates the extraction of signal features from multivariate time-series data stored in multiple Excel files. Key steps include:
- Segmentation of data for feature extraction.
- Extraction of signal features like spectral density, variability, and correlation metrics.
- Visualization of feature trends and comparative analysis.
- Classification of signals using extracted features.

## Features Extracted

The following statistical features are computed for each signal segment:
1. **SCCA**: Signal Cross-Correlation Analysis.
2. **SSVL**: Signal Sum of Variability Levels.
3. **TACR**: Temporal Amplitude Correlation Ratio.
4. **SDHC**: Sum of Data Horizontal Components.
5. **SH45 and SH135**: Signal Harmonics at 45° and 135°.
6. **SDTC**: Signal Distance-to-Centroid.
7. **SABP**: Sum of Absolute Bilateral Products.
8. **SCRA**: Signal Cross-Relation Amplitude.
9. **SHCA** and **SCTA**: Harmonic and Temporal Cross-Determinants.

These features capture spectral, amplitude, and spatial properties of the signal, making them useful for signal classification and analysis.

## Code Structure

The code is organized into the following sections:

- **Initialization and Data Loading**: 
  Sets paths and loads data from Excel files into MATLAB arrays.
  
- **Feature Extraction**: 
  The `ex_features` function calculates the statistical features for each data segment.
  
- **Visualization**: 
  Generates plots for each feature, allowing trend analysis across data segments.
  
- **Classification**: 
  Uses MATLAB’s `classificationLearner` for supervised learning on the extracted features.
  
- **Mean and Variance Analysis**: 
  Calculates and plots mean and variance of extracted features across different data sets.

## Requirements

- MATLAB (R2021a or newer recommended)
- Excel files containing the time-series data to be analyzed
- Signal Processing Toolbox (if using MATLAB’s signal processing functions)

## Usage

1. **Load the Data**: 
   Place the Excel files in the `dataPath` directory as specified in the code. Each file should contain multivariate time-series data, with features extracted from the second column.

2. **Run the Code**: 
   Execute the main script in MATLAB to perform feature extraction, analysis, and visualization.

3. **Output**: 
   Plots of extracted features will be saved in the `outputPath` directory, and statistical features will be displayed and compared in bar graphs. 

4. **Classification**: 
   Use the `classificationLearner` app to classify signals based on extracted features.

## Results and Visualizations

- **Feature Plots**: Time-series plots of each extracted feature show trends across different data segments.
- **Scatter Plots**: Compares signal data across columns.
- **Feature Comparison**: Bar graphs display comparative analysis of features from different data sets.
- **Mean and Variance Analysis**: Highlights variability within features, aiding in understanding signal consistency.

