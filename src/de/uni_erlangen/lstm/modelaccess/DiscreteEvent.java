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

import java.util.logging.Logger;

import org.apache.commons.math3.ode.events.EventHandler;


/**
 * G-Stop facility
 * ODE must be solved up to some target state is reached, 
 * with a known value of the state but an unknown occurrence time.
 * 
 * @author liampetti
 *
 */
public class DiscreteEvent implements EventHandler {
	public final static Logger LOGGER = Logger.getLogger(DiscreteEvent.class.getName());
	
	private int i;
	private double target;
	private double t;
	private boolean dirIncrease;
	
	public int getI() {
		return i;
	}

	public void setI(int i) {
		this.i = i;
	}

	public double getTarget() {
		return target;
	}

	public void setTarget(double target) {
		this.target = target;
	}

	public boolean isDirIncrease() {
		return dirIncrease;
	}

	public void setDirIncrease(boolean dirIncrease) {
		this.dirIncrease = dirIncrease;
	}

	/**
	 * Initialise this class with the required parameters
	 * 
	 * @param i 			Parameter
	 * @param target		The target value for the given parameter
	 * @param dirIncrease	Stop when increasing (false means stopping on decreasing value)
	 */
	public DiscreteEvent (int i, double target, boolean dirIncrease) {
		this.i = i;
		this.target = target;
		this.dirIncrease = dirIncrease;
	}
	
	public double getTime() {
		return t;
	}

	@Override
	public Action eventOccurred(double t, double[] y, boolean increasing) {
		if (increasing && !dirIncrease) {
			return EventHandler.Action.CONTINUE;
		} else {
			this.t = t;
			return EventHandler.Action.STOP;
		}
	}

	@Override
	public double g(double t, double[] y) {
		return y[i] - target;
	}

	@Override
	public void init(double arg0, double[] arg1, double arg2) {
		// TODO Auto-generated method stub
	}

	@Override
	public void resetState(double arg0, double[] arg1) {
		// TODO Auto-generated method stub
	}
}
