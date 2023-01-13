# SISREMO
Simple Storage Reqirement Model for a 100% Renewable Energy System

There are many researches dedicated to 100% renewable energy system modelling and their models are much more sophisticated than the model presented here. Some publications in this field are:

1. T. M. Clack etal., Evaluation of a proposal for reliable low-cost grid power with 100% wind, water, and solar, Proc Natl Acad Sci USA, 14, 16, (2017), www.pnas.org/cgi/doi/10.1073/pnas.1610381114
2.  Ziegler etal., Storage Requirements and Costs of Shaping Renewable Energy Toward Grid Decarbonization, Joule 3, 2134–2153 (2019), https://doi.org/10.1016/j.joule.2019.06.012
3.  Ch. Breyer, etal., On the History and Future of 100% Renewable Energy Systems Research, IEEE Access, 10, (2022), DOI: 10.1109/ACCESS.2022.3193402
4.  D. Bogdanov, Radical transformation pathway towards sustainable electricity via evolutionary steps, Nature Communications volume 10, Article number: 1077 (2019), https://www.nature.com/articles/s41467-019-08855-1
5.  O. Ruhnau, O. Qvist, Storage requirements in a 100% renewable electricity system: extreme events and inter-annual variability, Environ. Res. Lett. 17 044018, (2022), DOI 10.1088/1748-9326/ac4dc8

The model described here is intended to show the relationships between energy production, energy consumption and storage requirements to clarify a little. It is easier to understand than the models described in the literature and anyone who is interested can play with the code themselves with a little [Julia](https://julialang.org/) programming experience.

## Data

The data used are from the site [Energy-Charts API](https://api.energy-charts.info/)

1. download power data as json from https://api.energy-charts.info/power using the REST API:
    * execute: **load_ise_energy_chart_data(start_year, end_year)**
        in **storage/data_energy_charts.jl**
        the minimum start_year is 2015
2. Parse downloaded json files and store date in a **hdf5** file
   * execute **run_ise_json_to_hdf5(false, 2015, 2022)**
        in **storage/data_energy_charts.jl**

Energy-Chart data used in this model are:

1. Dates: UNIX timestamps are converted to Julia DateTime objects
2. Load
3. Sum of Wind offshore, Wind onshore and Solar

Load and renewable time series data
![RP](figures/RP.png)

The real data are adapted to mimic a 100% renewable eneryg system by detrending and scaling.

## Detrend and Scale Date

Detrending is done by fitting teh data to a second order polynomial
### Detrend Load

$k$  - polynomial order
$n$  - number of data points

1. Load trend, polynomial fit
$L_{t} = \operatorname{polynomial\_fit}(L, k)$
2. Detrend
$L_{d} = \dfrac{L_{t}[n/2]}{L_{t}} \; L$

$L$  - Load [MW]
$L_{t}$ -  trend of Load, poynomial fit
$L_{d}$ - detrended Load

Detrended Load data

![Load_d](figures/Load_detrended.png)

### Detrend and Scale Renewables

**function scale_and_detrend(Load::Vector{Float64}, RP::Vector{Float64})**

1. First scaling
$R_{s} = R \; \dfrac{\operatorname{mean}(L)}{\operatorname{mean}(R)}$

2. Renewables trend, polynomial fit
$R_t = \operatorname{polynomial\_fit}(R_s, k)$

3. Detrend
$R_d = \dfrac{R_t[n/2]}{R_t} R_s$

4. Scale again
$R_{ds} = R_d \dfrac{\operatorname{mean}(L_d)}{\operatorname{mean}(R_d)}$

5. Diffeernce between Renewables and Laod
$\Delta P = (R_{ds} - L_d)$

$R$ - renewable energy [MW]
$R_{s}$ -scaled $R$
$R_d$ - detrended $R$
$R_t$ - detrended and scaled $R$

Detrended and scaled renewable power data
![RP_d](figures/RP_detrended.png)

Differences beteeen reneable power and load
![RP_d](figures/RP_diff_detrended.png)

## Compute Storage Fill Level as Function of Time

Given storage capacity and an overproduction capacity factor the storage fille level is computed;

**compute_storage_level(dates, Load, RP, eunit, over_production, storage_capacity)**

The algorithm in short is:

$\Delta P = R - L$

if $\Delta P > 0$ and $S < S_{capacity}$

$\quad S = S + \Delta P \Delta t$

elseif $\Delta P < 0$ and $S > 0$

$\quad S = S - \Delta P \Delta t$

end

Storage fill level over time for different combinations of **storage capacity** and renewable overproduction factor **op**

![storage](figures/storage_fill.png)


