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

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;
import java.util.logging.Logger;

/**
 * Writes specified data to a CSV file
 * 
 * @author liampetti
 *
 */
public class CSVWriter {
	public final static Logger LOGGER = Logger.getLogger(CSVWriter.class.getName());
	
	/*
	 * Print one line of data
	 */
	public void WriteString(String filename, String output, boolean append) {				
		try {
			PrintStream fileStream = new PrintStream(new FileOutputStream(filename, append));
			fileStream.println(output);
			fileStream.close();
		} catch (IOException e) {
			LOGGER.severe(e.toString());
		} 
	}
	
	/*
	 * Print one line of data
	 */
	public void WriteArray(String filename, double[] outputs, boolean append) {				
		try {
			PrintStream fileStream = new PrintStream(new FileOutputStream(filename, append));
			String printer = "";
			for (int i=0;i<outputs.length;i++) {
				printer += outputs[i] + ";";
			}
			fileStream.println(printer);
			fileStream.close();
		} catch (IOException e) {
			LOGGER.severe(e.toString());
		} 
	}
	
	/*
	 * Print all data from a list
	 */
	public void WriteList(String filename, List<double[]> dataset, boolean append) {				
		try {
			// Set the fileoutput stream to append mode
			PrintStream fileStream = new PrintStream(new FileOutputStream(filename, append));			
			String printer = "";
					
			for (double[] data : dataset) {
				printer = data[0] + ";";
				for (int i=1;i<data.length;i++) {
					printer += data[i] + ";";
				}
				fileStream.println(printer);
			}
									
			fileStream.close();
		} catch (IOException e) {
			LOGGER.severe(e.toString());
		} 
	}
	
	/*
	 * Clear file
	 */
	public void Clear(String filename) {				
		try {
			PrintStream fileStream = new PrintStream(new FileOutputStream(filename, false));
			fileStream.print("");
			fileStream.close();
		} catch (IOException e) {
			LOGGER.severe(e.toString());
		} 
	}
}
