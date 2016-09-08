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

package de.uni_erlangen.lstm.file;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Logger;

/**
 * Reads a CSV data file for processing
 * 
 * @author liampetti
 *
 */
public class CSVReader {
	public final static Logger LOGGER = Logger.getLogger(CSVReader.class.getName());
	
	// File data
	private String filename;
	private BufferedReader br = null;
	private String splitter;
	
	private boolean finished = false;
	
	// Default settings
	public CSVReader() {
		filename = "input.csv";
		splitter = ";";
		initBuffer();
		
	}
	
	public CSVReader(String filename, String splitter) {
		this.filename = filename;
		this.splitter = splitter;
		initBuffer();
	}
	
	public void initBuffer() {
		try {
			br = new BufferedReader(new FileReader(filename));
		} catch (FileNotFoundException e) {
			LOGGER.severe(e.toString());
		}
	}
	
	public boolean finished() {
		return finished;
	}
	
	public String[] getStrings() {
		String[] currentList = new String[0];
		String line = "";
		
		// Load our input list from the CSV file (reload)
		try { 
			if ((line = br.readLine()) != null) {	 
				currentList = line.split(splitter); 
			} else {
				finished = true;
				br.close();
			}
		} catch (IOException e) {
			LOGGER.severe(e.toString());
		}
		
		return currentList;
	}

}
