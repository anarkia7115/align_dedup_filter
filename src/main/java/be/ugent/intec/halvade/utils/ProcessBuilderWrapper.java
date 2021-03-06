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

import be.ugent.intec.halvade.tools.ProcessException;
import java.io.*;
import java.util.Arrays;

/**
 *
 * @author ddecap
 */
public class ProcessBuilderWrapper {
    private class ProcMon implements Runnable {
        private final Process _proc;
        private int _complete;
        
        public ProcMon(Process proc) {
            _proc = proc; 
            _complete = 1;
        }
        
        public int getState() { return _complete; }

        @Override
		public void run() {
            try {
                int val = _proc.waitFor();
                Logger.DEBUG("process ended with " + val);
            } catch (InterruptedException ex) {
                _complete = -1;
            }
            _complete = 0;
        }
    }
    // the process itself
    Process p;
    String[] command;
    String libdir;
    // streams
    StreamGobbler stderr;
    StreamGobbler stdout;
    BufferedWriter stdin;
    FileOutputStream outStream = null;
    ProcMon mon;
    long startTime, estimatedTime;
    int threads = 1;
    
    public ProcessBuilderWrapper(String[] command, String libdir) {
        this.command = command;
        this.libdir = libdir;
    }
    
    public void setThreads(int threads) {
        this.threads = threads;
    }
    
    public void startProcess() throws InterruptedException {
        startProcess(null, null);
    }
    
    public void startProcess(boolean redirectStreams) throws InterruptedException {
        if(redirectStreams) {
            startProcess(System.out, System.err);
        } else {
            startProcess(null, null);
        }
    }
    
    public void startProcess(String redirectFile) throws InterruptedException, IOException {
    	this.outStream = new FileOutputStream(redirectFile, true);
        startProcess(outStream, System.err, true);
    }
    
    public void startProcess(FileOutputStream filestream, PrintStream stderr_, boolean redirectToFile) throws InterruptedException {
        try {
            Logger.DEBUG("running command " + Arrays.toString(command));            
            ProcessBuilder builder = new ProcessBuilder(command);
            if(libdir != null) {
                builder.environment().put(
                        "LD_LIBRARY_PATH",
                        "$LD_LIBRARY_PATH:" + libdir);
                builder.environment().put(
                         "PYTHONPATH",
                         "$PYTHONPATH:" + libdir);
            }
            builder.environment().put("CILK_NWORKERS", "" + threads);
            p = builder.start();
            mon = new ProcMon(p);
            Thread t = new Thread(mon);
            t.start();
            startTime = System.currentTimeMillis();
            if (filestream != null) {
            	// p.getInputStream returns the process's output stream !
            	this.stdout = new StreamGobbler(p.getInputStream(), filestream);
                this.stdout.start();
            }
            if(stderr_ != null) {
                this.stderr = new StreamGobbler(p.getErrorStream(), stderr_, "[PROCESS_ERR]");
                this.stderr.start();
            } 
            stdin = new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
        } catch (IOException ex) {
            Logger.EXCEPTION(ex);
            throw new ProcessException(ex.getMessage(), -1);
        }
        if(mon.getState() != 1) 
            Logger.DEBUG("command not started!: " + Arrays.toString(command));             
    }
    
    public void startProcess(PrintStream stdout_, PrintStream stderr_) throws InterruptedException {
        try {
            Logger.DEBUG("running command " + Arrays.toString(command));            
            ProcessBuilder builder = new ProcessBuilder(command);
            if(libdir != null) {
                builder.environment().put(
                        "LD_LIBRARY_PATH",
                        "$LD_LIBRARY_PATH:" + libdir);
                builder.environment().put(
                         "PYTHONPATH",
                         "$PYTHONPATH:" + libdir);
            }
            builder.environment().put("CILK_NWORKERS", "" + threads);
            p = builder.start();
            mon = new ProcMon(p);
            Thread t = new Thread(mon);
            t.start();
            startTime = System.currentTimeMillis();
            if (stdout_ != null) {
            	this.stdout = new StreamGobbler(p.getInputStream(), stdout_);
                this.stdout.start();
            }
            if(stderr_ != null) {
                this.stderr = new StreamGobbler(p.getErrorStream(), stderr_, "[PROCESS_ERR]");
                this.stderr.start();
            } 
            stdin = new BufferedWriter(new OutputStreamWriter(p.getOutputStream()));
        } catch (IOException ex) {
            Logger.EXCEPTION(ex);
            throw new ProcessException(ex.getMessage(), -1);
        }
        if(mon.getState() != 1) 
            Logger.DEBUG("command not started!: " + Arrays.toString(command));             
    }
    
    /**
     * 
     * @return 1 is running, 0 is completed, -1 is error 
     */
    public int getState() {
        return mon.getState();
    }
    
    public long getExecutionTime() {
        return this.estimatedTime;
    }
    
    public int getExitState() {
        return p.exitValue();
    }
    
    public boolean isAlive() {
        return mon.getState() >= 0;
    }
    
    public InputStream getSTDOUTStream() {
        return p.getInputStream();
    }
    
    public InputStream getSTDERRStream() {
        return p.getErrorStream();
    }
    
    public OutputStream getSTDINStream() {
        return p.getOutputStream();
    }
    
    public BufferedWriter getSTDINWriter() {
        return stdin;
    }
    
    public int waitForCompletion() throws InterruptedException, IOException {
        int val = p.waitFor();
        //if (outStream != null) outStream.close();
        Boolean outStreamNotNull = false;
        if (outStream != null) outStreamNotNull = true;
        Logger.DEBUG("outStream not null : " + outStreamNotNull.toString() );

        estimatedTime = System.currentTimeMillis() - startTime;
        return val;
    }
}
