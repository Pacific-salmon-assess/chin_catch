# Estimates of Chinook salmon condition and abundance from SWVI troll survey

Data and code associated with estimates of Chinook salmon condition and spatiotemporal distribution from DFO fisheries-independet troll survey on southwest Vancouver Island (2019-2023). Manuscript (**Seasonal variability in condition and spatial distribution of Chinook Salmon: Implications for ecosystem-based management**) in press at **Transactions of the American Fisheries Society**.

Scripts:
1. catch_sdm_mvrw.R: spatiotemporal models for different ecological groups
2. clean_age_data.R: clean aging data for exploratory figures
3. data_figs.R: data cleaning, composition and individual condition figures, individual condition GAMs
4. env_data.R: local sea surface temperature data from lighthouses and La Perouse buoy (used to evaluate seasonal/interannual variability)
5. maps.R: generate study area and raw catch maps
6. prep_prediction_grid.R: generate grid for spatial predictions
7. sampling_sim.R: sensitivity analysis evaluating alternative sampling design
8. theme_sleek2.R: helper function for ggplot theme
