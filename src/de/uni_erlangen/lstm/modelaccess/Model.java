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

package de.uni_erlangen.lstm.modelaccess;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.AdamsBashforthIntegrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import de.uni_erlangen.lstm.file.CSVWriter;
import de.uni_erlangen.lstm.models.adm1.DAEModel;
import de.uni_erlangen.lstm.models.adm1.DigesterParameters;
import de.uni_erlangen.lstm.models.adm1.StateVariables;

/**
 * Class for controlling the ADM1 model, can be run on a separate thread
 * 
 * @author liampetti
 *
 */
public class Model implements Runnable {
	public final static Logger LOGGER = Logger.getLogger(Model.class.getName());
	
	private double[] x;
	private double[] u;
	private double[] param;
	private double S_H_ion;
	private double start;
	private double end; 
	private List<DiscreteEvent> events;
	private boolean finished;
	private boolean onlineRecord; // Record model to CSV
	private double resolution; // How often to sample data from continuous model
	private double progress;
	private boolean dae;
	private double fix_pH;
		
	/**
	 * Initialise model using custom parameters and outputs
	 * 
	 * @param start 		Initial time 
	 * @param end 			Final time
	 * @param parameters 	Custom digester parameters
	 * @param outputs 		Custom outputs
	 * @param contModel 	Build a continuous output model of the solution
	 */
	public Model(double start, double end, DigesterParameters parameters, StateVariables initial, StateVariables influent, boolean onlineRecord) {
		dae = true;
		fix_pH = -1.0;
		this.onlineRecord = onlineRecord;
		this.resolution = 0.01041666667; // 15 minutes in days as standard resolution
		u = influent.getVar(); // Influent
		x = initial.getVar(); // Output (initial reactor conditions)
		param = parameters.getParameters();	
		init(start, end);
	}

	public void init(double start, double end) {
		this.events = new ArrayList<DiscreteEvent>();
		this.start = start;
		this.end = end;
		this.progress = start;
		
		// Flow rate is set by influent
		x[35] = u[35]; // Effluent flow rate = Influent flow rate
		
		// Initialise the S_H_ion
		double factor = (1.0/param[0] - 1.0/param[1])/(100.0*0.083145);
		double K_w = Math.pow(10,-param[2])*Math.exp(55900.0*factor); // T adjustment for K_w 
		double phi = x[24]+(x[10]-x[31])-x[30]-(x[29]/64.0)-(x[28]/112.0)-(x[27]/160.0)-(x[26]/208.0)-x[25];
		S_H_ion = (-phi*0.5)+0.5*Math.sqrt(phi*phi+(4.0*K_w)); // SH+	
	}
	
	public void setTime(double start, double end) {
		this.start = start;
		this.end = end;
		this.progress = start;
	}
	
	public void setInfluent(StateVariables influent) {		
		u = influent.getVar(); // Influent
		x[35] = u[35]; // Effluent flow rate = Influent flow rate
	}
	
	public void setInitial(StateVariables initial) {
		x = initial.getVar(); // Initial effluent
	}
	
	public void setParameters(DigesterParameters parameters) {		
		param = parameters.getParameters();	
	}
	
	public void setDAE (boolean dae) {
		this.dae = dae;
	}
	
	public void setpH (double ph) {
		this.fix_pH = ph;
	}
	
	public void addEvent (DiscreteEvent event) {
		this.events.add(event);
	}
	
	public void addEvents (List<DiscreteEvent> events) {
		this.events = events;
	}
	
	public void setResolution(double res) {
		this.resolution = res;
	}
	
	/**
	 * Run the model using set parameters
	 */
	public void simulate() {		
		finished = false;		
		/*
		 * Integrator selection 
		 */
		//FirstOrderIntegrator integrator = new HighamHall54Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		//FirstOrderIntegrator integrator = new DormandPrince54Integrator(1.0e-12, 100.0, 1.0e-12, 1.0e-12);
		//FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		//FirstOrderIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		FirstOrderIntegrator integrator = new AdamsBashforthIntegrator(2, 1.0e-14, 100.0, 1.0e-10, 1.0e-10);
		//FirstOrderIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		
		// influent values, digester parameters, S_H_ion, dae system
		final DAEModel ode = new DAEModel(u, param, S_H_ion, dae, fix_pH);
		//FirstOrderDifferentialEquations ode = model; 
		
		// Records progress
		StepHandler progHandler = new StepHandler() {
		    public void init(double t0, double[] y0, double t) {
		    }
		            
		    public void handleStep(StepInterpolator interpolator, boolean isLast) {
		    	progress = interpolator.getCurrentTime();		        
		    }
		};
		integrator.addStepHandler(progHandler);
		
		/*
		 * Continuous model recorded in CSV
		 */
		if (onlineRecord) {
			final CSVWriter writer = new CSVWriter();
			StepHandler stepHandler = new StepHandler() {
				double prevT = 0.0;
				
			    public void init(double t0, double[] y0, double t) {
			    }
			            
			    public void handleStep(StepInterpolator interpolator, boolean isLast) {
			        double   t = interpolator.getCurrentTime();
			        if (t-prevT > resolution) {
			        	// Add time to the beginning of the array
						double[] timemodel = new double[ode.getDimensions().length+1];
						timemodel[0] = t;
						
						// We need to pull variables (S_h2 and acid-base) directly from the model if using DAE
						for (int i=1;i<timemodel.length;i++) {
							timemodel[i] = ode.getDimensions()[i-1];
						}
							
			        	writer.WriteArray("cont_model_output.csv", timemodel, true);
			        	prevT = t;
			        }
			    }
			};
			integrator.addStepHandler(stepHandler);
		}
		
		/*
		 * Add event handlers for discrete events
		 * maxCheck - maximal time interval between switching function checks (this interval prevents missing sign changes in case the integration steps becomes very large)
    	 * conv - convergence threshold in the event time search
    	 * maxIt - upper limit of the iteration count in the event time search
		 */
		if (events.size() > 0) {
			for (DiscreteEvent event : events) {
				double maxCheck = Double.POSITIVE_INFINITY;
				double conv = 1.0e-20;
				int maxIt = 100;
				integrator.addEventHandler(event, maxCheck, conv, maxIt);
			}
		}
			
		integrator.integrate(ode, start, x, end, x);

		/*
		 * Return the time that the discrete event occurred
		 */
		if (events.size() > 0) {
			for (DiscreteEvent event : events) {
				if (event.getTime() < end) {
					end = event.getTime();
				}
			}
		}

		// We need to pull variables (S_h2 and acid-base) directly from the model
		x = ode.getDimensions();
		
		finished = true;
	}
	
	public boolean isFinished() {
		return finished;
	}
	
	public double getProgress() {
		return progress;
	}
	
	public double[] getX() {
		return x;
	}
	
	public void setX(double[] x) {
		this.x = x; // Initial effluent
	}
	
	public double[] getU() {
		return u;
	}
	
	public double[] getParam() {
		return param;
	}
	
	public double getEnd() {
		return end;
	}
	
	/**
	 * Allows the simulation to run on a separate thread
	 * 
	 * @see java.lang.Runnable#run()
	 */
	@Override
	public void run() {		
		simulate();	
	}
}
