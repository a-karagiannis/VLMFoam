# VLMFoam

Hello and thanks for visiting this corner of GitHub :-)

The C++ code in this directory constitutes a custom numerical solver for the open-source, finite-volume CFD platform OpenFOAM, which I compiled to simulate relatively simple cases of pure water vapour condensation in conditions of supersonic expansion. Examples of such conditions would be the flow in a supersonic nozzle or over the blades of a steam turbine. I used it in the less common context of miniature nozzles for water-based spacecraft propulsion systems. These are called Vapourising Liquid Microthrusters (VLM), hence the solver's name. If interested in more details, you could find them in my master thesis [1].

I validated the code against the experiments of Moore et al (1973) [2], a widely used benchmark for low-pressure condensation studies of this type. I have included one of these validation cases as an example in this repository. Keep in mind that the modelling of metastable steam condensation is notorious for the uncertainty involved, so the solver's performance may deteriorate when compared to other experiments. A good example of such inconsistencies can be seen for other codes in literature, in [3].

As usual in OpenFOAM, you can use the solver by simply copying the files to your local working directory and compiling them. I have checked that the solver compiles successully in OpenFOAM 6, on my personal laptop's 18.04 Ubuntu Linux. Perhaps adaptations may be needed for earlier/later versions, as portions of the underlying structure of OpenFOAM tend to change between releases, sometimes significantly. 

#### A word of caution
This code was my first step with C++ programming and OpenFOAM and I have not had time to review it since I first made it. As a result, I advise caution for the possibility of poor programming practices, perhaps even mistakes. I welcome any constructive feedback/criticism. I myself was originally planning to expand the code such that it can simulate mixtures of two gases rather than just H2O vapour, as well as simulate condensation of substances other than just water, such as hydrocarbons. Unfortunately, a change in personal circumstances now makes it unlikely I will be able to do any of this before early 2022, if at all. It is also unlikely I will be able to process any potential pull requests until then.

Be that as it may, I will be happy if this solver turns out to be of use to someone other than just myself, so I welcome you to use it and modify it as you wish. 

Cheers! 

## The solver

The solver is an extension of the already existing "rhoCentralFoam" solver for supersonic flows, included by default in the standard OpenFOAM release (as of 2020). Details about the underlying numerics of "rhoCentralFoam" are outlined in [4]. Details about the code in this repository can be found in my master thesis, in [1].


## Custom thermodynamic model

While the thermophysics of H2O steam at low pressures is not dramatically different from an ideal gas's, the modelling of condensation requires a real-gas model, as well as models for the properties of the resulting liquid and the properties at saturation. Normally, an OpenFOAM solver should be a stand-alone entity and the thermodynamics should be implemented separately into OpenFOAM's internal themodynamics framework. I found it more convenient to implement the models for the liquid and saturation properties directly in the solver itself and implement only the real-gas model in OpenFOAM's internal structure. You can find the real-gas model in the "src" folder and the rest in the solver's directory.

During my thesis I also implemented a simple table-look-up method for the thermodynamics of high-pressure steam. If I find some time, I will clean it up and add it here.

Regardless of whether you use the included thermodynamic package or one of your own construction, you should keep in mind to check:
l. that your chosen EoS and the rest of thermodynamic relations are thermodynamically consistent with each other and
l. that your EOS exhibits thermodynamic stability at the region of the phase diagram that is of interest to your case.


## Custom boundary conditions

Due to some particularities in the flow cases I was treating, I have also implemented a couple of custom boundary conditions (BCs), which you can find in the respective folder in the solver. Unless you are treating similar flow physics, these boundary conditions are probably irrelevant to your work and you may ignore this section.

I mentioned earlier that I wrote this solver because I needed to check for condensation inside micronozzles. The gaseous flow in these minute devices exhibits (at least) three complications from a modelling perspective:
l. Extremely rapid expansions that cause instabilities in the numerics of "rhoCentralFoam"
l. A portion of the flow domain downstream of the nozzle outlet has to be included for the flow solution to be realistic [5]. This means that "rhoCentralFoam", an explicit solver, has to compute a flow that goes from subsonic in the nozzle convergent to supersonic in the nozzle divergent and back to subsonic at the far outlet. Since the solver has not really been designed for that, it poses problems in the realisticity of the solution and the solution speed. 
l. A small degree of rarefaction in the flow, especially near the nozzle exit.

To partially alleviate the first two issues, I implemented a custom "subsonic-supersonic pressure outlet" boundary condition for the pressure. This was adapted from the work of Kraposhin [6]. To partially alleviate the last issue, I coded slip boundary conditions at the walls for the velocity and temperature, after the work of Gokcen [7]. 


## References

[1] [Master thesis](https://repository.tudelft.nl/islandora/object/uuid%3Aabe34357-627e-4df4-92eb-73e590ab79a6?collection=education)

[2] Moore, M. J., Walters, P. T., Crane, R. I., and Davidson, B. J. (1973). Predicting the fog-drop size inwet-steam turbines. InProceedings of the Conference on Wet Steam 4, Coventry, UK. Institution ofMechanical Engineers. Paper C37/73

[3] Starzmann et al. (2018). Results of the international wet steam modeling project. Proceedings of the Institution of Mechanical Engineers, Part A: Journal of Power and Energy,232(5):550–570. DOI: 10.1177/0957650918758779

[4] Greenshields, C. J., Weller, H. G., Gasparini, L., and Reese, J. M. (2009). Implementation of semi-discrete, non-staggered central schemes in a colocated, polyhedral, finite volume framework, forhigh-speed viscous flows.International Journal for Numerical Methods in Fluids 

[5] Ivanov, M., Markelov, G., Ketsdever, A., and Wadsworth, D. (1999). Numerical study of cold gas mi-cronozzle flows. In 37th Aerospace Sciences Meeting and Exhibit. American Institute of Aeronautics and Astronautics. DOI: 10.2514/6.1999-166

[6] Kraposhin, M., Bovtrikova, A., and Strijhak, S. (2015).   Adaptation of Kurganov-Tadmor nu-merical scheme for applying in combination with the PISO method in numerical simulation of flows in a wide range of mach numbers.Procedia Computer Science, 66:43–52.  DOI:10.1016/j.procs.2015.11.007

[7] Gokcen, T., Maccormack, R., and Chapman, D. (1987). Computational fluid dynamics near the con-tinuum limit. In8th Computational Fluid Dynamics Conference. American Institute of Aeronauticsand Astronautics. DOI: 10.2514/6.1987-1115
