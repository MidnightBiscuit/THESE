# fluo-variations_optimal-temp

I had a code working to compute P, dP/dT, and find the temperature for which fluo was maximum. In the version of Dec 2, commit 100e6ce, the code provides a satisfying result.

## Computation of Probabilities

$\texttt{MB}$ is the Maxwell-Boltzmann distribution. It provides the probability to find an atom with a given velocity in a gas with temperature T.

$\texttt{pfl\_dop} = \rho_{ee}$ is the atomic ray profile with Doppler effect. It provides the probability of excitation given Rabi frequency, detuning, lambda and velocity of atom. It is considered equal to the excited population, i.e the proportion of excited atoms.

The fluorescence F is related to the product of both probabilities $ F \propto \texttt{MB}\times \texttt{pfl\_dop} = \texttt{prob\_fluo}$. This product is still a probability, lets call it $P$. $dP/dT = 0$ gives us the temperature for which the fluo is maximum, $d^2P/dT^2 = 0$ gives us the temperature for which the fluorescence variation  $dP/dT$ is maximum.

The idea is to search for the setting where the variation of fluorescence is the highest for a given input in energy. The energy introduced by the GiantMolecule is brought under the form of thermal energy. Thus we want to set the cloud in a configuration with the best variation of fluorescence for a variation in temperature. There are two ways of achieving this. First considering infinitesimal differences, this requires to compute $dP/dT$. Second considering macroscopic differences $\Delta P$ and $\Delta T$. Because we know how much energy the GiantMolecule transfers to the cloud, we can compute the expected $\Delta T$. Then search for the optimal $\Delta P$ given the determines $\Delta T$. In the simulation article, the GMol with incident energy 50eV, transfers 50meV, which produces a 100mK increase of temperature in the cloud.

## Code version from Dec 2 100e6ce

The Doppler limit temperature depending the detuning, everything is ok with this result.

<img src="/home/adrian/Documents/Programmes/Python/THESE/Notes/Cop_Lim_Temp.png" alt="Cop_Lim_Temp" style="zoom:50%;" />

The extrema from $dE/dT$ are retrieved. Here it appears OK. In the simulation $\delta = \Gamma$ and the temperature for maximum fluo is $\approx 500$ mK, which is given by this graph.

<img src="/home/adrian/Documents/Programmes/Python/THESE/Notes/Bifurcation-1.png" alt="Bifurcation-1" style="zoom:50%;" />

<img src="/home/adrian/Documents/Programmes/Python/THESE/Notes/Ideal_temperature_fluo_max.png" alt="Ideal_temperature_fluo_max" style="zoom: 50%;" />

The extremum are retrieved from data showed below.

<img src="/home/adrian/Documents/Programmes/Python/THESE/Notes/Integral_derivate_+_extrema.png" alt="Integral_derivate_+_extrema" style="zoom:50%;" />

## Code from Dec 08

It appears i forgot a square root in MB function, line 9. The right result is `return (m_Ca/(np.pi*2*kb*T))**(1/2) * np.exp(-m_Ca*v**2/(2*kb*T)) ` But i had written `**(1)`.

It is okay now with the square root.



## Some computations compared to data from article

In article conditions, I determine :

- the temperature for which the fluorescence is maximum. I compute the produce $P = \texttt{MB}\times\texttt{pfl\_dop}$ and see when $\mathrm{d}P/\mathrm{d}T = 0$ or at least when $\mathrm{d}P/\mathrm{d}T$ is minimum. I found in the data that fluo is max at 168 mK, but the computation gives me 225 mK. Maybe it is due to the fact that during the GMol is passing, the MB is not very accurately describing the gas ?
- the temperature for which the fluo goes back to its original level $P_{cold}$. I try to find when $P = P_{cold}$ after the maximum. In the data it is for 3-4 K. In the computation i cannot find a satisfying value.
- the fluo variation is maximum when $\mathrm{d}^2P/\mathrm{d}T^2 = 0$