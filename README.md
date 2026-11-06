# Differential Responses to a Bright and Noisy World: How Sensory Pollution Influences Morphology and Condition Across Three Songbird Species

**Authors**: 

## Abstract

1.	As anthropogenic land use expands and fundamentally restructures the global biosphere, the pervasive infiltration of artificial light at night and human-generated noise introduces novel selective pressures that penetrate even the most isolated habitats. These sensory pollutants rarely operate in isolation, and their combined influence on adult avian phenotypes remains a poorly understood dimension of global environmental change.
2.	Using a landscape-scale experimental system, we investigated how independent and interactive effects of macro- and microhabitat light and noise influence adult morphology and condition across a community of three cavity-nesting songbirds.
3.	Our results reveal that sensory pollutants drive divergent phenotypic outcomes across species rather than uniform shifts. Exposure to sensory pollution generated size- and age-based spatial sorting, conditional penalties, sex-specific sensitivities, and, under combined exposures, synergistic conditional benefits.
4.	These species spanned a distinct functional continuum of responses to sensory pollution. Mountain Bluebirds (Sialia currucoides) incurred noise-driven, size-based spatial sorting alongside independent mass and condition penalties in both macro- and microhabitat light and noise. In contrast, Western Bluebirds (Sialia mexicana) demonstrated apparent conditional benefits, with older adults concentrated in co-polluted habitats with synergistic body condition gains under microhabitat dual-pollutant exposure. Finally, Ash-throated Flycatchers (Myiarchus cinerascens) displayed distinct sensory asymmetry, showing acoustic resilience while incurring microhabitat light-driven condition penalties across both sexes, wherein longer-winged females occupied the brightest habitats despite suffering steeper condition declines.
5.	Broadly, our findings demonstrate that the physical toll of sensory change – manifested through species-specific responses in adult morphology, age structure, and body condition – depends on spatial scale and pollutant combinations, highlighting the need to account for individual sensory sensitivities and life-history constraints to predict avian persistence in an increasingly bright and noisy world.


## Repository Structure:
All R code and data within the repository can be found in the corresponding folders.

- `Data`: Contains all data used in the manuscript.
  - `Adult_Morphology_Cavity.csv`: Individual-level data on adult morphology. Contains columns for species (MOBL – Mountain Bluebird, WEBL – Western Bluebird, ATFL – Ash-throated Flycatcher), biological sex, individual service band ID, color band ID, sensory pollution treatment (Control, Light, Noise, Light+Noise), cluster section, site, date of collection, year, nest ID, age (SY – Second Year, AHY – After Hatch Year, ASY – After Second Year), reproductive stage (TE – territory establishment, Eggs – has eggs, Chx – has nestlings), lux (lux light values), LAeq (A-weighted noise levels), mass, tarsus length, wing length, individual measurer initials.

- `Figures`: Contains all R code to replicate figures in the manuscript.
    - `Manuscript_Figures.csv`: Script to reproduce Figures 1 - 3 found in the Main Text and Figure S1 of the Supplementary Materials using `Adult_Morphology_Cavity.csv`.
      
- `R Code`: Contains all R scripts of main analyses in the manuscript:
   - `ATFL_Adult_Morph.R`: Uses `Adult_Morphology_Cavity.csv` to analyze Ash-throated Flycatcher morphology. Runs species-specific linear and linear mixed-effects models across four response variables (wing length, tarsus length, body mass, and body condition calculated as a scaled mass index). Performs a two-step model selection approach: first determining the parsimonious random effect structure via AICc, and second evaluating competing macrohabitat (treatment groups) and continuous microhabitat (lux and LAeq) fixed-effect models with interactions, controlling for covariates.
   - `MOBL_Adult_Morph.R`: Uses `Adult_Morphology_Cavity.csv` to analyze Mountain Bluebird morphology. Implements the same two-step model selection framework (random structure optimization followed by information-theoretic fixed-effect ranking) across wing length, tarsus length, body mass, and body condition under macro- and microhabitat sensory pollutant treatments.
   - `WEBL_Adult_Morph.R`: Uses `Adult_Morphology_Cavity.csv` to analyze Western Bluebird morphology. Evaluates wing length, tarsus length, body mass, and body condition across the two-step model selection framework while incorporating age (SY, AHY, ASY) as an additional fixed covariate to account for experience variance.
