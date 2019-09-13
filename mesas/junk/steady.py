# -*- coding: utf-8 -*-
"""Storage selection (SAS) functions: example with one flux out at steady state

Runs the rSAS old_model for a synthetic dataset with one flux in and out
and steady state flow

Theory is presented in:
Harman, C. J. (2015), Time-variable transit time distributions and transport:
Theory and application to storage-dependent transport of chloride in a watershed,
Water Resour. Res., 51, doi:10.1002/2014WR015707.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sas_functions import Piecewise
from sas_model import Model, SASTimeseries

from mesas import sas

plt.ion()
makeplots = True
# =====================================
# Generate the input timeseries
# =====================================
# length of the dataset
S_0 = 10.  # <-- volume of the uniformly sampled store
Q_0 = 1.0  # <-- steady-state flow rate
T_0 = S_0 / Q_0
N = 30
n_substeps = 5
# Steady-state flow in and out for N timesteps
data_df = pd.DataFrame(index=range(N))
data_df['J'] = np.ones(N) * Q_0
data_df['Q'] = np.ones(N) * Q_0
# A timeseries of concentrations
data_df['tracer'] = np.ones(N)
# =========================
# Parameters needed by SAS
# =========================
solute_parameters = {
    'tracer': {}
}
# =========================
# Create the SAS functions
# =========================
# Parameters for the rSAS function
# The uniform distribution extends between S_T=a and S_T=b.
sas_fun_Q = Piecewise(npiece=5, ST_max=S_0)
sas_ts = {
    'Q': SASTimeseries(sas_fun_Q, N=N),
}
# =============
# Run the old_model
# =============
# Run it
model = Model(
    data_df=data_df,
    sas_ts=sas_ts,
    solute_parameters=solute_parameters,
    n_substeps=3,
    verbose=True,
)
model.run()
outputs = model.result
# %%
# Timestep-averaged outflow concentration
# ROWS of C_Q are t - times
# COLUMNS of PQ are q - fluxes
C_Qm1 = outputs['C_Q'][:, 0, 0]
# Age-ranked storage
# ROWS of ST are T - ages
# COLUMNS of ST are t - times
# LAYERS of MS are s - solutes
ST = outputs['ST']
# Timestep-averaged backwards TTD
# ROWS of PQ are T - ages
# COLUMNS of PQ are t - times
# LAYERS of PQ are q - fluxes
PQm = outputs['PQ'][:, :, 0]
# Timestep-averaged outflow concentration
# ROWS of C_Q are t - times
# COLUMNS of PQ are q - fluxes
# Use SAS.transport to convolve the input concentration with the TTD
C_Qm2, C_mod_raw, observed_fraction = sas.transport(
    PQm, data_df['tracer'], model.solute_parameters['tracer']['C_old'])

if makeplots:
    # ==================================
    # Plot the age-ranked storage
    # ==================================
    print('Plotting ST at the last timestep')
    # The analytical solution for the age-ranked storage is
    T = np.arange(N+1)
    ST_exact = S_0 * (1 - np.exp(-T/T_0))
    # plot this with the SAS estimate
    fig = plt.figure(1)
    plt.clf()
    plt.plot(ST[:, -1], 'b-', label='SAS old_model', lw=2)
    plt.plot(ST_exact, 'r-.', label='analytical solution', lw=2)
    plt.ylim((0, S_0))
    plt.legend(loc=0)
    plt.ylabel('$S_T(T)$')
    plt.xlabel('age $T$')
    plt.title('Age-ranked storage')
    #
    # =====================================================================
    # Outflow concentration estimated using several different TTD
    # =====================================================================
    # Lets get the instantaneous value of the TTD at the end of each timestep
    print('Getting the instantaneous TTD')
    PQi = np.zeros((N+1, N+1))
    PQi[:, 0] = sas_ts['Q'](ST[:, 0], 0)
    PQi[:, 1:] = np.r_[[sas_ts['Q'](ST[:, i+1], i) for i in range(N)]].T
    # Lets also get the exact TTD
    print('Getting the exact solution')
    n = 100
    T = np.arange(N*n+1.)/n
    PQe = np.tile(1-np.exp(-T/T_0), (N*n+1, 1)).T
    # Use the transit time distribution and input timeseries to estimate
    # the output timeseries for the exact and instantaneous cases
    print('Getting the concentrations')
    C_Qi, C_mod_raw, observed_fraction = sas.transport(
        PQi, data_df['tracer'], model.solute_parameters['tracer']['C_old'])
    C_Qei, C_mod_raw, observed_fraction = sas.transport(
        PQe, data_df['tracer'].repeat(n), model.solute_parameters['tracer']['C_old'])
    # This calculates an exact timestep-averaged value
    C_Qem = np.reshape(C_Qei, (N, n)).mean(axis=1)
    # Plot the results
    print('Plotting concentrations')
    fig = plt.figure(2)
    plt.clf()
    plt.step(np.arange(N), C_Qem, 'r', linestyle='-',
             label='mean exact', lw=2, where='post')
    plt.step(np.arange(N), C_Qm1, 'g', linestyle='--',
             label='mean SAS internal', lw=2, where='post')
    plt.step(np.arange(N), C_Qm2, 'm', linestyle=':',
             label='mean SAS.transport', lw=2, where='post')
    plt.plot((np.arange(N*n) + 1.)/n, C_Qei, 'r-', label='inst. exact', lw=1)
    plt.plot(np.arange(N)+1, C_Qi, 'b:o', label='inst. SAS.transport', lw=1)
    plt.legend(loc=0)
    plt.ylim((0, 1))
    plt.ylabel('Concentration [-]')
    plt.xlabel('time')
    plt.title('Outflow concentration')
    plt.draw()
    plt.show()