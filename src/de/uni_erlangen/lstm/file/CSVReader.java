/*
 * jADM1 -- Java Implementation of Anaerobic Digestion Model No 1
 * ===============================================================
 *
 * Copyright 2016 Liam Pettigrew
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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
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
	private String splitter;
	private BufferedReader br = null;
	
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
	
	/**
	 * Gets the last line from the given CSV file
	 * 
	 * @return A string array of the last line
	 */	
	public String[] getLastStrings() {
		File file = new File(filename);
	    RandomAccessFile fileHandler = null;
	    try {
	    	fileHandler = new RandomAccessFile(file, "r" );
	        long fileLength = fileHandler.length() - 1;
	        StringBuilder sb = new StringBuilder();

	        for(long filePointer = fileLength; filePointer != -1; filePointer--){
	            fileHandler.seek(filePointer);
	            int readByte = fileHandler.readByte();

	            if(readByte == 0xA) {
	                if(filePointer == fileLength) {
	                    continue;
	                }
	                break;
	            } else if(readByte == 0xD) {
	                if(filePointer == fileLength-1) {
	                    continue;
	                }
	                break;
	            }
	            sb.append((char) readByte);
	        }

	        String lastLine = sb.reverse().toString();
	        String[] currentList = lastLine.split(splitter);
	        return currentList;
	    } catch(Exception e) {
	    	LOGGER.severe(e.toString());
	        return null;
	    } finally {
	        if (fileHandler != null )
	            try {
	                fileHandler.close();
	            } catch (IOException e) {
	            	LOGGER.warning(e.toString());
	            }
	    }
	}
	
	/**
	 * Gets the next line from the given CSV file
	 * 
	 * @return A string array of from the current line, starting at the beginning
	 */
	public String[] getNextString() {
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
