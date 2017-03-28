/*
 * jADM1 -- Java Implementation of Anaerobic Digestion Model No 1
 * ===============================================================
 *
 * Copyright 2015 Liam Pettigrew
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * ********************************************************************************************
 */

package de.uni_erlangen.lstm.main;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import de.uni_erlangen.lstm.file.CSVReader;
import de.uni_erlangen.lstm.file.CSVWriter;
import de.uni_erlangen.lstm.modelaccess.DiscreteEvent;
import de.uni_erlangen.lstm.modelaccess.Model;
import de.uni_erlangen.lstm.models.adm1.BSM2Defaults;
import de.uni_erlangen.lstm.models.adm1.DigesterParameters;
import de.uni_erlangen.lstm.models.adm1.StateVariables;

/**
 * Main class allows user access to the model through a command line interface
 * 
 * Command line arguments ->
 * -steady	Run steady state simulation
 * -dynamic Run dynamic simulation
 * -cont 	Write continuous output model to CSV file
 * -s  		Start time (in days)
 * -f		Finish time (in days)
 * -in		Influent filename for steady (one line) or dynamic (multiple lines)
 * -init	Reactor initial conditions filename
 * -param	Reactor parameters
 * -step 	Step size for dynamic model influent (in days)
 * -ode 	Run ODE model (very slow!)
 * -event 	Add state event to the simulation to tell it when to stop, three variables: variable number, variable value, rising/falling (true/false)
 * 
 * @author liampetti
 * 
 */
public class Main {
	public final static Logger LOGGER = Logger.getLogger(Main.class.getName());
	
	private Model model;
	private boolean steady;
	private String[] args;
	
	private double stime;
	private double start; // Model start time
	private double finish; // Model end time
	private StateVariables initial; // Digester initial conditions
	private StateVariables influent; // Influent
	private DigesterParameters parameters; // Model parameters
	private boolean modOut; // Store all model outputs (needed for plotting)
	private double step; // Adjust time step size for model outputs
	private boolean dae; // Tells the model to run the algebraic equations
	private List<DiscreteEvent> events; // Discrete event detection
	private CSVReader dynamicIn; // Input file for dynamic influent

	public void start(String[] args) {
		this.args = args;
		boolean spec = false; 
		if (args.length > 0) {
			for (int i=0;i<args.length;i++) {
				switch (args[i]) {
					case "-steady":		runSteady();
										steady = true;
										spec = true;
										break;
					case "-dynamic": 	runDynamic();
										steady = false;
										spec = true;
										break;
					default:			break;
				}
			}	
		} 
		
		// User did not specify what simulation type, default to steady
		if (!spec) {
			runSteady();
			steady = true;
		}
	}
	
	private void runSteady() {	
		CSVWriter writer = new CSVWriter();
		stime = System.currentTimeMillis();
		events = new ArrayList<DiscreteEvent>();
		// Setup model outputs and parameters (default is BSM2)
		BSM2Defaults defaults = new BSM2Defaults();
		initial = new StateVariables();
		initial.setVar(defaults.DigesterInit());
		influent = new StateVariables();
		influent.setVar(defaults.Influent());
		parameters = new DigesterParameters();
		// No command line arguments, run a default setup
		start = 0.0;
		finish = 200.0;
		step = 0.01041666667; // 15 minutes in days as standard resolution
		modOut = false;
		dae = true;
		
		checkArgs();

		model = new Model(start, finish, step, parameters, initial, influent, modOut, "steady_out.csv");	
		model.setDAE(dae);		
		model.addEvents(events);
		
		if (modOut) {
			writer.Clear("cont_model_output.csv");
		}
		
		new Thread(model).start();
		
		while(!model.isFinished()) {
			try {
				Thread.sleep(3000); // Wait 3 seconds before checking
				System.out.println("Progress = " +
						String.format("%.2f",(model.getProgress()/finish)*100)
						+ "%");
			} catch (InterruptedException e) {
				LOGGER.severe(e.toString());
			} 
			
		}
		
		double[] x = model.getX();
		double[] u = model.getU();
		
		String output = "Simulation time; " + (System.currentTimeMillis()-stime) + 
				"; Start; " + start + 
				"; Finish; " + model.getEnd() + "\n";
		for (int i=0;i<x.length;i++) {
			output += "State no; " + (i+1) + 
					";\t Influent; " + u[i] + 
					";\t Effluent; " + x[i] + "\n";
	 	}
		System.out.println(output);
		writer.WriteString("steady_result.csv", output, true);
	}
	
	private void runDynamic() {
		double stime = System.currentTimeMillis();
		CSVWriter writer = new CSVWriter();
		writer.Clear("dynamic_output.csv");
		events = new ArrayList<DiscreteEvent>();
		// Setup model outputs and parameters (default is BSM2)
		BSM2Defaults defaults = new BSM2Defaults();
		initial = new StateVariables();
		initial.setVar(defaults.DigesterInit());
		influent = new StateVariables();
		dynamicIn = new CSVReader("digesterin.csv", ",");
		parameters = new DigesterParameters();
		// No command line arguments, run a default setup
		start = 0.0;
		finish = 609.0;
		step = 0.01041666667; // 15 minutes in days	
		// Continuous output models are not used in dynamic models at the moment
		modOut = false;	
		dae = true;
		
		checkArgs();

		model = new Model(start, start+step, step, parameters, initial, influent, modOut, "dynamic_out.csv");
		model.setDAE(dae);
		model.addEvents(events);
		
		int t = 0;
		
		while (!dynamicIn.finished()) {
			String[] inString = dynamicIn.getNextString();		
			if (inString.length > 0) {
				double[] in = new double[inString.length];
				for (int i=0;i<in.length;i++) {
					in[i] = Double.parseDouble(inString[i]);
				}
				influent.setVar(in);
				model.setInfluent(influent);
			}
			
			model.setTime(start, start+step);
			model.run();
			
			// Add time to the beginning of the array and save to csv
			double[] timemodel = new double[model.getX().length+1];
			timemodel[0] = start;
			for (int i=1;i<timemodel.length;i++) {
				timemodel[i] = model.getX()[i-1];
			}
			writer.WriteArray("dynamic_output.csv", timemodel, true);
			
			start = start+step;
			if (t%(Math.round(finish/100)) == 0) {
				System.out.println("Progress = " + String.format("%.2f",(start/finish)*100) + "%");
			}
			t++;
		}
		
		System.out.println("Simulation time; " + (System.currentTimeMillis()-stime));
	}
	
	private void checkArgs() {
		if (args.length > 0) {
			for (int i=0;i<args.length;i++) {
				switch (args[i]) {
					case "-cont":	modOut = true;
									break;
					case "-s":		start = Double.parseDouble(args[i+1]);
									break;
					case "-f":		finish = Double.parseDouble(args[i+1]);
									break;
					case "-in":		if (steady) {	
										influent.readVar(args[i+1]);
									} else {
										dynamicIn = new CSVReader(args[i+1], ",");
									}
									break;
					case "-init":	initial.readVar(args[i+1]);
									break;
					case "-param": 	parameters.readParameters(args[i+1]);
									break;
					case "-step":	step = Double.parseDouble(args[i+1]);
									break;
					case "-ode":	dae = false;
									break;
					case "-event":	DiscreteEvent event = new DiscreteEvent(Integer.parseInt(args[i+1]),
										Double.parseDouble(args[i+2]),
										Boolean.parseBoolean(args[i+3]));
									events.add(event);
									break;
					default:		break;
				}
			}
		}
	}
	
	public static void main(String[] args) {
		Main main = new Main();
		main.start(args);		
	}
}
