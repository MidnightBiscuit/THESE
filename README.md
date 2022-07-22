# Codes for my PhD THESIS

Those are the codes I used during my thesis in order to compute parameters, analyse data from simulations and experiments, produce plots. In this document, I would like to highlight some of the programs according to their particular relevance. It is appropriate in the case the program deals with a very important feature of my thesis, and is clear and can serve as an example or for a discussion. A program can also be highlighted because it is related to a critical point or a new aspect that our team is not familiar with. The featured highlighted programs are indicated with a $\star$ so they are more clearly visible. Keep in mind that some programs I used may not be present in the list, although I did my best to gather and summarize everything I used. The references to sections, figures and tables are done according to my thesis.

[TOC]

## Laser cooling of trapped ions

### Laser cooling

$\star$ `fluo-variations_optimal-temp.ipynb` and `fluo-variations_optimal-temp_minimal` are dedicated to the study of the ion fluorescence when submitted to laser cooling. This program accounts for single ion and ion ensemble. See subsections 4.2. and 4.3 and figures 4.3 to 4.6. Inside the code many details are provided about the computations.

```Rabi.ipynb``` computes Rabi frequency and saturation intensity for a given condition.

```Saturation_ions_par_laser.ipynb``` is related to the fluorescence saturation measurements. See subsection 4.7.1.

```scan_397_fluo.ipynb``` and ```scan_866_fluo.ipynb``` are codes dedicated to the analysis of the ion cloud fluorescence measured during a scan of a laser (397 or 866). See subsection 4.7.2.

### Ion trapping

#### General purpose

$\star$ ```2021_Mathieu_parameters.ipynb``` computes all trapping parameters for GiantMol and TADOTI, given trapping voltages.

```axial_potential.ipynb``` computes the axial trapping voltage. Quadratic and gaussian. See subsection 1.4.2, figure 1.6.

$\star$ ```Continued fraction.ipynb``` computes the $\beta_u(a_u,q_u)$ for a given $q_u$ in both adiabatic and non-adiabatic cases (with the continued fraction). Also computes for a given $q$ several other parameters such as $U_{RF}$, $f_x$ ... Also computes the axial voltage required to have a unity aspect ratio $\alpha = 1$.

<<<<<<< HEAD
=======
$\star$ ```computation_V_f_q_piege_nonadia.ipynb``` computes the trapping frequencies and Mathieu parameters as a function of the radio-frequency field amplitude. It considers non-adiabatic case and adiabatic case. It provides a beautiful and informative graph.

>>>>>>> abcec4537c18105df704f4925240f6224d0cf0c5
`Potentiels_piege.ipynb` computes trapping potential in 1D, 2D, for quadrupole and other multipoles.

```trap_basic_equations.ipynb``` is a small code to compute trapping parameters (frequencies, aspect ratio, maximum axial potential for stable trapping) from trapping voltages. See subsection 1.1.2 and figures 1.2 and 1.3.

#### Image analysis

$\star$ `Compensation_Image_analysis.ipynb` . This programs do many things. First it is devoted to the signal processing of camera images according to the method discussed subsection 3.4.2. It is also dedicated to the contact potential determination. First this program loads camera images (from Tucsen), apply several processing transformations such as

- gaussian spatial filter
- threshold filtering and mask creation
- edge detection
- ellipse fitting (see Box 3.4.1) (see around # Plot the least squares ellipse in this program)

See figure 2.7 and 3.2. This program is also used in order to measure contact potential. The center of the ellipse is retrieved and a plot of its position as a function of trapping parameter can be done. Ultimately, the axial contact potential value is provided. See subsection 2.4.1 and figure 2.6.

#### Tickle analysis

<<<<<<< HEAD
`data_tickle.ipynb` and `data_tickle_reverse.ipynb` retrieve fluorescence measured during a tickle experiment and plot the fluo. Those programs also provide a peak detection in the fluorescence in order to automatically detect the main tickle frequency. The program associates the right frequency with the peak, provided the required elements are provided (scan amplitude, trapping voltages and so on). See subsection 2.3.1 and figure 2.2.
=======
`data_tickle.ipynb` and `data_tickle_reverse.ipynb` retrieve fluorescence measured during a tickle experiment and plot the fluo. Those programs also provide a peak detection in the fluorescence in order to automatically detect the main tickle frequency. The program associates the right frequency with the peak, provided the required elements are provided (scan amplitude, trapping voltages and so on). See subsection 2.3.1 and figure 2.2. See also ```data_tickle-220713.ipynb``` and ```data_tickle-220713.ipynb```
>>>>>>> abcec4537c18105df704f4925240f6224d0cf0c5

`data_UDC.ipynb` and `data_URF.ipynb` are also related to tickle measurement. The purpose of those programs is to study the secular frequency as a function of $U_{DC}$ or $U_{RF}$. See subsection 2.3.2 and figures 2.3, 2.4.

## The molecular source

### Hydrodynamics

$\star$ ```Capillary_flow.ipynb``` computes capillary flow according to several models : Bernoulli, Hagen-Poiseuille and Wutz/Adam model. See Wissdorf article for the later model. See subsection 5.2.2.

```Pipe_flow.ipynb``` computes flow parameters in the capillary according to two models pointed out by Wissdorf.

### Electrostatics and charged particle optics

```Analytic potential.ipynb``` computes the analytic potential from a source fed with an electrostatic potential. Not very useful.

```Bender_curve.ipynb``` computations about the particle trajectories in electrostatic quadrupolar bender. See subsection 6.2.6, equations 6.22 and 6.23.

$\star$ `Charged_particle_optics_matrix.ipynb` matrix computation for electrostatic lens properties. Provides the focal length depending on the electrode potential. See subsections 6.2.1 to 6.2.4 and figure 6.1.

`Electrodes_potential.ipynb` computes electric field generated by point charges. Same remark as for `Analytic potential.ipynb`.

```Potential_CSV_Agros2D.ipynb``` retrieves potentials saved with Agros2D (as .csv files) and plot them. See subsection 6.2.4, figure 6.3.

#### MCP analysis

`MCP_ES_plot-190404` and `MCP_ES_plot-191105` are used in order to process data measured by oscilloscope in MCP. See subsection 7.2.

## Numerical simulations

$\star$ `Langevin-init-Quench-simple-ARTICLE-OPEN-ALL_STATS.ipynb` is one of the numerous program dedicated to the analysis of numerical simulations. This program is processing the data from the simulation of the collision of a GiantMol with a trapped and laser cooled ion cloud. It is adapted to the version of the program starting with a Langevin initialisation. It handles the statistical apsects. Good luck because this program is truly humongous. See sections 9.1, 9.2 and 9.3.

### RF_Temp_Fit

Codes to work the RF Heating simulations 20211201, 20220218.

`fit_T_Good.ipynb` small code to analyze data from RF Heating simulations. This is in fact the first code used to fit analytically the temperature curves. It is outdated. You would rather use code in the directory `RF_Temp_fit`, such as the following.

$\star$ `20211201_RF_Heating.ipynb` is dedicated to the analysis of RF Heating simulations. The date is corresponding to the date to which the data are labelled in the simulation server.  See the whole section 9.4.

## Other

`Egun_carac.ipynb`processes electron current measurements carried out in GiantMol. See subsection 4.5.2 and table 4.2, figure 4.12. Also this program integrates considerations about the black body radiation and an analytical attempt to directly correlate the electron current to the filament current. The maximum current the filament can admit is estimated. See Box 4.5.1.

`lock_carac.ipynb` is intended to process laser intensities measured with photodiodes in and out the doubling cavity. See subsection 4.4.3 and figures 4.9 and 4.10.

`MB-distribution.ipynb` just works out Maxwell-Boltzmann distribution and produces some figures, along with numerical values for given conditions..
