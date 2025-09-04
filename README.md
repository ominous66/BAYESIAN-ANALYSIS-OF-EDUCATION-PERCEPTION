# Bayesian-Analysis-of-Education-Perception | Jan’25–Apr’25**

## Project Overview

This project applies Bayesian quantile analysis to examine public perceptions of the educational value of online versus in-person courses. Drawing on a nationally representative U.S. dataset, the study aims to identify demographic and socioeconomic factors that shape attitudes toward online education, providing nuanced insights for educators and policymakers.

## Motivation

The rapid expansion of online education—accelerated by the COVID-19 pandemic—has raised important questions about its perceived value compared to traditional classroom teaching. While online learning offers flexibility and broader access, public opinion remains divided, especially among different demographic groups. Understanding these perceptions is crucial for designing effective, inclusive educational policies and programs.

## Data and Methodology

- **Dataset:**  
  - Sourced from the Pew Social Trends and Demographics Project (2011), with 24 variables and 1,591 cleaned observations after excluding missing and ambiguous responses.
  - Variables include age, gender, education, employment status, region, race, income, prior online course experience, and more.

- **Statistical Approach:**  
  - Developed a **Bayesian Markov Chain Monte Carlo (MCMC)** algorithm for binary quantile regression, using 12,500 chains and 2,500 burn-in iterations.
  - Modeled binary outcomes (perceived equal value of online vs. in-person courses) via latent responses, simulated through truncated normal distributions.
  - Quantile regression framework captures heterogeneity in preferences across the latent utility scale, overcoming limitations of traditional mean-based models like probit or logit.

## Key Findings

- **Preference for Online Education:**
  - **Older individuals:** +2% higher probability of valuing online education.
  - **Full-time workers:** +8.6% higher probability.
  - **Prior online experience:** +10.9% higher probability.
  - **Females:** Show greater propensity for online education, especially at lower and middle quantiles[1][2][3].

- **Lower Preference:**
  - **Highly educated individuals (Post-bachelor’s):** −15.1% lower probability of valuing online education compared to those with a high school education or less.
  - **Bachelor’s degree holders:** Also show a negative effect, though less pronounced.

- **Other Insights:**
  - **Regional variation:** Notable differences, with the Northeast and South less likely to value online courses compared to the Midwest.
  - **No significant effects:** Race, income, and area of residence (urban/suburban/rural) do not significantly influence perceptions, suggesting broad accessibility and a narrowing digital divide.

## Policy Implications

- Online education is unlikely to fully replace traditional classrooms, especially for highly educated populations who value in-person interaction.
- Online programs should be tailored to the needs of older adults, full-time professionals, and those with prior online experience, leveraging the flexibility these formats offer.
- The democratization of access—evidenced by the lack of race and income effects—highlights online education’s potential to reach diverse populations, but quality and rigor must be maintained.

## Technical Skills & Tools

- **Statistical Modeling:** Bayesian quantile regression, MCMC, latent variable modeling.
- **Programming:** R (bqror package), Python, data wrangling and visualization.
- **Data Analysis:** Pandas, NumPy, statistical inference, and result interpretation.

## References

- Ojha, M., & Rahman, M. A. (2021). *Do online courses provide an equal educational value compared to in-person classroom teaching? Evidence from U.S. survey data using quantile regression*. Education Policy Analysis Archives.

## Acknowledgments

Project completed as part of ECO545 under the guidance of **Prof. Arshad Rahman** at the Indian Institute of Technology Kanpur.
