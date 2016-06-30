/*
 * Copyright (C) 2014 ddecap
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package be.ugent.intec.halvade.utils;

import java.io.*;

/**
 *
 * @author ddecap
 */
public class StreamGobbler extends Thread {
    InputStream is;
    PrintStream stream;
    FileOutputStream filestream;
    String prefix;

    public StreamGobbler(InputStream is, PrintStream stream) {
        this.is = is;
        this.stream = stream;
        this.filestream = null;
        this.prefix = "";
    }
    
    public StreamGobbler(InputStream is, PrintStream stream, String prefix) {
        this.is = is;
        this.stream = stream;
        this.filestream = null;
        this.prefix = prefix;
    }

    public StreamGobbler(InputStream inputStream, FileOutputStream filestream) {
        this.is = inputStream;
        this.stream = null;
        this.filestream = filestream;
        this.prefix = "";
	}

	@Override
    public void run() {
        try {
            if (filestream != null) {
            	final byte[] buffer = new byte[64 * 1024];//buffer size is: 64k
                int n = 0;
                while (-1 != (n = is.read(buffer))) {
                	filestream.write(buffer, 0, n);
                }
//            	DataInputStream isr = new DataInputStream(is);
//            	int index;
//            	while ((index = isr.read()) != -1) {
//            		filestream.write(index);
//            	}
            	filestream.flush();
            } else {
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);
            	String line = null;
                while ((line = br.readLine()) != null)
                    stream.println(prefix + "[" + Timer.getGlobalTime() + "] " + line);
            }
           
        }
        catch (IOException ioe) {
            Logger.EXCEPTION(ioe);
        }
    }
}
