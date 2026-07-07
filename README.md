Analysis Overview
The project includes four R scripts that perform the following operations:
- Data loading from CSV files using relative paths.
- Estimation of Alpha and Beta coefficients using Ordinary Least Squares (OLS) regression.
- Evaluation of non-linear relationships using Kernel smoothing regression.
- Diagnostic tests on OLS residuals, including density plots, Q-Q plots, Shapiro-Wilk, Jarque-Bera, and D'Agostino tests.
- Hypothesis testing using the Neyman-Pearson framework and Fisher tests for coefficient significance.
- Evaluation of test robustness through Statistical Power Analysis (Power Curves).
- Specific hypothesis testing to determine whether the stocks are cyclical or non-cyclical (testing Beta against a threshold of 1).

Project Structure
- `data/`: Contains the historical price CSV files.
- `scripts/`: Contains the R scripts for ENI and SHELL across both time periods.

How to Run
1. Download or clone this repository.
2. Open any script from the `scripts/` folder in RStudio.
3. Set your Working Directory to the script location (Session -> Set Working Directory -> To Source File Location).
4. Run the script to automatically load the data and generate the plots.
