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

package be.ugent.intec.halvade.tools;

import be.ugent.intec.halvade.hadoop.mapreduce.HalvadeCounters;
import java.io.IOException;
import java.io.InputStream;
import org.apache.hadoop.mapreduce.Mapper;
import be.ugent.intec.halvade.utils.*;
import java.net.URISyntaxException;
import java.util.ArrayList;

import org.apache.hadoop.mapreduce.Mapper.Context;

/**
 *
 * @author ddecap
 */
public class BWAMemInstance extends AlignerInstance {

	private static BWAMemInstance instance;
	private ProcessBuilderWrapper pbw;
	private SAMStreamHandler ssh;
	private ArrayList<String> as = new ArrayList<String>();
	private int retryTime = 0;
	private final int MAX_RETRY = 1;

	/**
	 * 
	 * This BWA instance runs BWA from stdin (custom provided BWA is needed)
	 */
	private BWAMemInstance(Context context, String bin) throws IOException, URISyntaxException {
		super(context, bin);
		String taskid = context.getTaskAttemptID().toString();
		taskid = taskid.substring(taskid.indexOf("m_"));
		ref = HalvadeFileUtils.downloadBWAIndex(context, taskid);
	}
	
	public int getAsSize() {
		return as.size();
	} 

	public int feedLine(String line) throws IOException {
		// gjx save line to mem
		as.add(line);
		return feedLine(line, pbw);
	}

	public void loadMapContextInMem(ProcessBuilderWrapper proc) {
		// if (proc.getState() != 1) {
		// Logger.DEBUG("writing \'" + line +"\' to process with state " +
		// proc.getState());
		// throw new IOException("Error when writing to process with current
		// state " + proc.getState());
		// }
		if (as.size() <=0){
			Logger.DEBUG("as size: " + Integer.toString(as.size()));
			return;
		}
		try {
			for (String line : as) {

				proc.getSTDINWriter().write(line, 0, line.length());

				proc.getSTDINWriter().newLine();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	protected void startAligner(Mapper.Context context) throws IOException, InterruptedException {
		// make command
		String customArgs = HalvadeConf.getCustomArgs(context.getConfiguration(), "bwa", "mem");
		String[] command = CommandGenerator.bwaMem(bin, ref, null, null, isPaired, true, threads, customArgs);
		pbw = new ProcessBuilderWrapper(command, bin);
		// run command
		// needs to be streamed to output otherwise the process blocks ...
		pbw.startProcess(null, System.err);
		// check if alive.
		if (!pbw.isAlive())
			throw new ProcessException("BWA mem", pbw.getExitState());
		pbw.getSTDINWriter();
		// make a SAMstream handler
		ssh = new SAMStreamHandler(instance, context, false);
		ssh.start();
	}

	/**
	 * 
	 * @return 1 is running, 0 is completed, -1 is error
	 */
	@Override
	public int getState() {
		return pbw.getState();
	}

	@Override
	public void closeAligner() throws InterruptedException, IOException {
		try {
			// close the input stream
			pbw.getSTDINWriter().flush();
			pbw.getSTDINWriter().close();
			// gjx
			Logger.DEBUG("***flush finished");
		} catch (IOException ex) {
			// gjx
			Logger.DEBUG("***in catch ex");
			Logger.EXCEPTION(ex);
		}

		int error = pbw.waitForCompletion();
		context.getCounter(HalvadeCounters.TIME_BWA_MEM).increment(pbw.getExecutionTime());
		if (error != 0) {
			// gjx
			Logger.DEBUG("***in err not 0");
			// if (retryTime < MAX_RETRY) {
			// Logger.DEBUG("retry " + Integer.toString(retryTime) + " times");
			// retryTime++;
			// startAligner(context);
			// closeAligner();
			// }
			// else {
			// throw new ProcessException("BWA mem", error);
			// }
			throw new ProcessException("BWA mem", error);
		}
		if (retryTime < MAX_RETRY) {
			Logger.DEBUG("retry " + Integer.toString(retryTime) + " times");
			retryTime++;
			startAligner(context);
			loadMapContextInMem(pbw);
			closeAligner();
		}
		if (retryTime != 0) {
			// in recursive
			return;
		}

		ssh.join();
		instance = null;
	}

	static public BWAMemInstance getBWAInstance(Mapper.Context context, String bin)
			throws IOException, InterruptedException, URISyntaxException {
		if (instance == null) {
			instance = new BWAMemInstance(context, bin);
			instance.startAligner(context);
		}
		AlignerInstance.context = context;
		return instance;
	}

	@Override
	public InputStream getSTDOUTStream() {
		return pbw.getSTDOUTStream();
	}

	@Override
	public void flushStream() {
		try {
			// close the input stream
			pbw.getSTDINWriter().flush();
		} catch (IOException ex) {
			Logger.EXCEPTION(ex);
		}
	}
}
