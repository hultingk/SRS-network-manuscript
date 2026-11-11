# SRS plant-pollinator network manuscript

### Scripts
**00-functions.R**: Script containing custom functions used through analysis.\
**01-sampling-completeness.R**: Script for estimating sampling completeness of plant-pollinator networks. \
**02-prepare-webs.R**: Script for wrangling networks into correct format for analysis.\
**03-hpc-bootstrap.R**: Script for calculating network metrics using bootstrapnet. \
**03b-network-d.R**: Script for calculating d'.\
**03c-diversity-metrics.R**: Script for calculating diversity metrics.\
**03d-SI-diversity-plots.R**: Script for analysis and plotting of supplementary diversity analysis.\
**04-network-metrics-wrangle.R**: Script for wrangling product of 03-hpc-bootstrap.R into correct format for analysis. \
**05a-diversity-glmmm.R**: Script for analyzing diversity metrics. \
**05b-bootstrap-network-glmm.R**: Script for analyzing network metrics. \
**06-interaction-turnover.R**: Script for calculating species turnover and interaction rewiring.

### Data
#### L1_data
**SRS-plant-pollinator.csv**
Each row contains a single observation of a plant-pollinator interaction. 
| Variable          | Description                                                                                                                                                                             |
| :---------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| block   | Site name (8, 10, 52, 57, 53N, 53S, or 54S) |
| patch         | Patch Type Identifier (B = Connected, W = Unconnected)    |
| location              | Grid location of plant-pollinator interaction within patch    |
| date              | Date of observation   |
| sampling_round              | One of seven sampling rounds |
| start_time         | Start time of plant-pollinator survey            |
| am_pm         | AM or PM of start time            |
| temp         | Air temperature in C at start of survey            |
| wind         | Categorical description of wind speed at start of survey            |
| cloud         | Estimated percent cloud cover at start of survey            |
| flower_species         | Identity of flower species that pollinator was observed interacting with   |
| pollinator_species         | Identity of pollinator species observed   |
| order         | Order of pollinator species observed  |
| family         | Family of pollinator species observed  |
| genus         | Genus of pollinator species observed  |
| species         | Specific epithet of pollinator species observed  |

#### L2_boot_metrics
**h2.RData**: H2' metric bootstrapping result\
**links.RData**: Average number of links per species bootstrapping result\
**net_d_mean.csv**: d' of higher and lower level result\
**niche.RData**: Niche overlap of higher and lower level bootstrapping result\
**nodf.RData**: Nestedness bootstrapping result\
**robustness.RData**: Robustness to secondary extinction bootstrapping result\
**shannon.RData**: Shannon diversity of interactions bootstrapping result\
**webs.matrix.RData**: networks for analysis

#### L3_metrics
**diversity_metrics.csv**: result of diversity metric calculations \
**network_metrics_boot.csv**: result of network metric calculations












