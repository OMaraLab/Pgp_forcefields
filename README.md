# Pgp_forcefields

Materials used in simulating and analysing the simulations of P-glycoprotein in a POPC/CHOL bilayer in 5 different force fields: AMBER 99SB-ILDN, CHARMM 36, OPLS-AA/L, GROMOS 54A7, MARTINI 2.2

Materials are available in data/.

Analysis code used is available in analysis/.

## Parameter origins

The AMBER 99SB-ILDN force field was published in [1]. Lipid parameters were obtained from SLipids, at [2].


The CHARMM 36 force field was published in [3]. The lipid membrane was built with CHARMM-GUI ([4], [5]) using lipid parameters from [6].

The OPLS-AA/L force field was published in [7]. Lipid parameters were obtained from [8] and [9].

The GROMOS 54A7 force field was published in [10]. Lipid parameters were published in [11] and [12].

The MARTINI 2.2 force field was published in [13]. The bilayer was generated with insane [14].

## Analysis

Analysis was conducted with cpptraj [15] and Python code. The Python code was built on MDAnalysis [16][17] and MDTraj [18]. The seaborn [19] and Plotly [20] libraries were used in visualisation. 

All material is provided under an MIT license. Please cite the above references if you use this material in published work.

## References

[1]: Lindorff-Larsen, K.; Piana, S.; Palmo, K.; Maragakis, P.; Klepeis, J. L.; Dror, R. O.; Shaw, D. E. *Improved Side-Chain Torsion Potentials for the Amber Ff99SB Protein Force Field.* Proteins 2010, 78 (8), 1950–1958. [https://doi.org/10.1002/prot.22711](https://doi.org/10.1002/prot.22711).

[2]: Jämbeck, J. P. M.; Lyubartsev, A. P. *Derivation and Systematic Validation of a Refined All-Atom Force Field for Phosphatidylcholine Lipids.* J. Phys. Chem. B 2012, 116 (10), 3164–3179. [https://doi.org/10.1021/jp212503e.](https://doi.org/10.1021/jp212503e)

[3]: Best, R. B.; Zhu, X.; Shim, J.; Lopes, P. E. M.; Mittal, J.; Feig, M.; MacKerell, A. D. *Optimization of the Additive CHARMM All-Atom Protein Force Field Targeting Improved Sampling of the Backbone ϕ, ψ and Side-Chain Χ1 and Χ2 Dihedral Angles.* J. Chem. Theory Comput. 2012, 8 (9), 3257–3273. [https://doi.org/10.1021/ct300400x.](https://doi.org/10.1021/ct300400x)

[4]: Jo, S.; Cheng, X.; Lee, J.; Kim, S.; Park, S.-J.; Patel, D. S.; Beaven, A. H.; Lee, K. I.; Rui, H.; Park, S.; Lee, H. S.; Roux, B.; MacKerell, A. D.; Klauda, J. B.; Qi, Y.; Im, W. *CHARMM-GUI 10 Years for Biomolecular Modeling and Simulation: Biomolecular Modeling and Simulation.* J. Comput. Chem. 2017, 38 (15), 1114–1124. [https://doi.org/10.1002/jcc.24660.](https://doi.org/10.1002/jcc.24660)

[5]: Jo, S.; Kim, T.; Iyer, V. G.; Im, W. CHARMM-GUI: *A Web-Based Graphical User Interface for CHARMM.* J. Comput. Chem. 2008, 29 (11), 1859–1865. [https://doi.org/10.1002/jcc.20945.](https://doi.org/10.1002/jcc.20945)

[6]: Klauda, J. B.; Venable, R. M.; Freites, J. A.; O’Connor, J. W.; Tobias, D. J.; Mondragon-Ramirez, C.; Vorobyov, I.; MacKerell, A. D.; Pastor, R. W. Update of the CHARMM All-Atom Additive Force Field for Lipids: Validation on Six Lipid Types. J. Phys. Chem. B 2010, 114 (23), 7830–7843. [https://doi.org/10.1021/jp101759q.](https://doi.org/10.1021/jp101759q)

[7]: Kaminski, G. A.; Friesner, R. A.; Tirado-Rives, J.; Jorgensen, W. L. *Evaluation and Reparametrization of the OPLS-AA Force Field for Proteins via Comparison with Accurate Quantum Chemical Calculations on Peptides.* J. Phys. Chem. B 2001, 105 (28), 6474–6487. [https://doi.org/10.1021/jp003919d.](https://doi.org/10.1021/jp003919d)

[8]: Maciejewski, A.; Pasenkiewicz-Gierula, M.; Cramariuc, O.; Vattulainen, I.; Rog, T. *Refined OPLS All-Atom Force Field for Saturated Phosphatidylcholine Bilayers at Full Hydration.* J. Phys. Chem. B 2014, 118 (17), 4571–4581. [https://doi.org/10.1021/jp5016627.](https://doi.org/10.1021/jp5016627)


[9]: Kulig, W.; Pasenkiewicz-Gierula, M.; Róg, T. *Topologies, Structures and Parameter Files for Lipid Simulations in GROMACS with the OPLS-Aa Force Field: DPPC, POPC, DOPC, PEPC, and Cholesterol.* Data in Brief 2015, 5, 333–336. [https://doi.org/10.1016/j.dib.2015.09.013.](https://doi.org/10.1016/j.dib.2015.09.013)

[10]: Schmid, N.; Eichenberger, A. P.; Choutko, A.; Riniker, S.; Winger, M.; Mark, A. E.; van Gunsteren, W. F. *Definition and Testing of the GROMOS Force-Field Versions 54A7 and 54B7.* Eur Biophys J 2011, 40 (7), 843–856. [https://doi.org/10.1007/s00249-011-0700-9.](https://doi.org/10.1007/s00249-011-0700-9)

[11]: O’Mara, M. L.; Mark, A. E. *The Effect of Environment on the Structure of a Membrane Protein: P-Glycoprotein under Physiological Conditions.* J. Chem. Theory Comput. 2012, 8 (10), 3964–3976. [https://doi.org/10.1021/ct300254y.](https://doi.org/10.1021/ct300254y)

[12]: Poger, D.; Van Gunsteren, W. F.; Mark, A. E. *A New Force Field for Simulating Phosphatidylcholine Bilayers.* J. Comput. Chem. 2010, 31 (6), 1117–1125. [https://doi.org/10.1002/jcc.21396.](https://doi.org/10.1002/jcc.21396)

[13]: de Jong, D. H.; Singh, G.; Bennett, W. F. D.; Arnarez, C.; Wassenaar, T. A.; Schäfer, L. V.; Periole, X.; Tieleman, D. P.; Marrink, S. J. *Improved Parameters for the Martini Coarse-Grained Protein Force Field.* J. Chem. Theory Comput. 2013, 9 (1), 687–697. [https://doi.org/10.1021/ct300646g](https://doi.org/10.1021/ct300646g).

[14]: Wassenaar, T. A.; Ingólfsson, H. I.; Böckmann, R. A.; Tieleman, D. P.; Marrink, S. J. *Computational Lipidomics with Insane: A Versatile Tool for Generating Custom Membranes for Molecular Simulations.* J. Chem. Theory Comput. 2015, 11 (5), 2144–2155. [https://doi.org/10.1021/acs.jctc.5b00209](https://doi.org/10.1021/acs.jctc.5b00209).

[15]: Roe, D. R.; Cheatham, T. E. *PTRAJ and CPPTRAJ: Software for Processing and Analysis of Molecular Dynamics Trajectory Data.* J Chem Theory Comput 2013, 9 (7), 3084–3095. [https://doi.org/10.1021/ct400341p](https://doi.org/10.1021/ct400341p).

[16]: Gowers, R. J.; Linke, M.; Barnoud, J.; Reddy, T. J. E.; Melo, M. N.; Seyler, S. L.; Domański, J.; Dotson, D. L.; Buchoux, S.; Kenney, I. M.; Beckstein, O. *MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics Simulations*. Proceedings of the 15th Python in Science Conference 2016, 98–105. [https://doi.org/10.25080/Majora-629e541a-00e](https://doi.org/10.25080/Majora-629e541a-00e).

[17]: Michaud-Agrawal, N.; Denning, E. J.; Woolf, T. B.; Beckstein, O. *MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.* J. Comput. Chem. 2011, 32 (10), 2319–2327. [https://doi.org/10.1002/jcc.21787](https://doi.org/10.1002/jcc.21787).

[18]: McGibbon, R. T.; Beauchamp, K. A.; Harrigan, M. P.; Klein, C.; Swails, J. M.; Hernández, C. X.; Schwantes, C. R.; Wang, L.-P.; Lane, T. J.; Pande, V. S. *MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories.* Biophysical Journal 2015, 109 (8), 1528–1532. [https://doi.org/10.1016/j.bpj.2015.08.015](https://doi.org/10.1016/j.bpj.2015.08.015).

[19]: Waskom, M. L. *Seaborn: Statistical Data Visualization.* Journal of Open Source Software 2021, 6 (60), 3021. [https://doi.org/10.21105/joss.03021.](https://doi.org/10.21105/joss.03021)

[20]: Inc, P. T. *Collaborative data science* [https://plot.ly.](https://plot.ly)
