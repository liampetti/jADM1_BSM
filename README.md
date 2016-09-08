#jADM1
A Java implementation of the Anaerobic Digestion Model No 1 (ADM1).

Based on the Matlab code produced for the Benchmark Simulation Model No 2 (BSM2) produced by the International Water Association (IWA) Task Group on Benchmarking of Control Strategies for Wastewater Treatment Plants. The original ADM1 was developed by the IWA Task Group for Mathematical Modelling of Anaerobic Digestion Processes.


This work has been created at:

Lehrstuhl für Strömungsmechanik,  
Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU)


Usage of this model can be found in the following conference papers:

* Pettigrew, L., Hubert, S., Groß, F., Delgado, A. (2015). Implementation of Dynamic Biological Process Models into a Reference Net Simulation Environment. In 16. ASIM Dedicated Conference Simulation in Production and Logistics. Dortmund: Fraunhofer IRB Verlag.
* Pettigrew, L., Delgado, A. (2016). Neural Network Based Reinforcement Learning Control for Increased Methane Production in an Anaerobic Digestion System. In 3rd IWA Specialized International Conference „Ecotechnologies for Wastewater Treatment“ (ecoSTP16). Cambridge, UK.


#Interface and Usage
A command line interface can be used to run simulations with the following command line arguments:

* -steady		  
  * Run steady state simulation, default setting uses parameters from BSM2 implementation
* -dynamic 			
  * Run dynamic simulation, default setting requires "digesterin.csv" exported from the BSM2 simulation (609 days with 15 minute intervals)
* -cont "filename" 	
  * Write continuous output model to CSV file
* -s 0.0			 	
  * Start time (in days)
* -f 0.0				
  * Finish time (in days)
* -in "filename"		
  * Influent filename for steady (one line) or dynamic (multiple lines)
* -init "filename"		
  * Reactor initial conditions filename
* -param "filename"	
  * Reactor parameters filename
* -step 0.0			
  * Step size for dynamic model influent (in days)
* -ode 				
  * Run as ODE (very slow!)
* -event 0 0.0 true 	
  * Add state event to the simulation to tell it when to stop, three variables: variable number, variable value, rising/falling (true/false)
  

For example, the default BSM2 200-day ADM1 steady state simulation can be run using the command 

> 				java -jar jADM1.jar -steady
 
 
 
CSV files contain the variable values in a ';' separated string. The influent file can contain 42 or 26 variables (BSM2 export).

0) monosaccharides (kg COD/m3)

1) amino acids (kg COD/m3)

2) total long chain fatty acids (LCFA) (kg COD/m3)

3) total valerate (kg COD/m3)

4) total butyrate (kg COD/m3)

5) total propionate (kg COD/m3)

6) total acetate (kg COD/m3)

7) hydrogen gas (kg COD/m3)

8) methane gas (kg COD/m3)

9) inorganic carbon (kmole C/m3)

10) inorganic nitrogen (kmole N/m3)

11) soluble inerts (kg COD/m3)

12) composites (kg COD/m3)

13) carbohydrates (kg COD/m3)

14) proteins (kg COD/m3)

15) lipids (kg COD/m3)	
	
16) sugar degraders (kg COD/m3)

17) amino acid degraders (kg COD/m3)

18) LCFA degraders (kg COD/m3)

19) valerate and butyrate degraders (kg COD/m3)

20) propionate degraders (kg COD/m3)

21) acetate degraders (kg COD/m3)

22) hydrogen degraders (kg COD/m3)

23) particulate inerts (kg COD/m3)


>24) cations (metallic ions, strong base) (kmole/m3)

>25) anions (metallic ions, strong acid) (kmole/m3)

>26) is actually Sva- = valerate (kg COD/m3)

>27) is actually Sbu- = butyrate (kg COD/m3)

>28) is actually Spro- = propionate (kg COD/m3)

>29) is actually Sac- = acetate (kg COD/m3)

>30) bicarbonate (kmole C/m3)

>31) ammonia (kmole N/m3)

>32) hydrogen concentration in gas phase (kg COD/m3)

>33) methane concentration in gas phase (kg COD/m3)

>34) carbon dioxide concentration in gas phase (kmole C/m3)


24/35) flow rate (m3/d) - Set by influent

25/36) temperature (deg C) - Set by digester

>37) methane volume (m3/d)

>38) gas volume (m3/d)

>39) pH

>40) optional

>41) optional


#Dependencies
* Requires the Apache Commons Mathematics Library 3.5


#Credits

Special thanks to Ulf Jeppsson, Christian Rosen and Darko Vrecko for use of their Matlab code of the ADM1, developed when (around 2004) they were all working together at the Department of Industrial Electrical Engineering and Automation (IEA), Lund University, Sweden.


Referenced works:

* Batstone, D.J., Keller, J., Angelidaki, I., Kalyuzhnyi, S.V., Pavlostathis, S.G., Rozzi, A., Sanders, W.T.M., Siegrist, H. & Vavilin, V.A. (2002). Anaerobic Digestion Model No. 1 (ADM1). IWA Scientific and Technical Report No 13, IWA Publishing, London, UK.
* Gernaey, K.V., Jeppsson, U., Vanrolleghem, P.A. & Copp, J.B. (2014). Benchmarking of Control Strategies for Wastewater Treatment Plants. IWA Scientific and Technical Report No. 23, IWA Publishing, London, UK.


#License
>Copyright 2016 Liam Pettigrew
>
> Licensed under the Apache License, Version 2.0 (the "License");
> you may not use this file except in compliance with the License.
> You may obtain a copy of the License at


>			http://www.apache.org/licenses/LICENSE-2.0        


> Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
> See the License for the specific language governing permissions and limitations under the License.

