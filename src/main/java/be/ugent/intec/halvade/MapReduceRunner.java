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

package be.ugent.intec.halvade;

import org.seqdoop.hadoop_bam.SAMRecordWritable;
import be.ugent.intec.halvade.hadoop.datatypes.ChromosomeRegion;
import be.ugent.intec.halvade.hadoop.mapreduce.HalvadeTextInputFormat;
import be.ugent.intec.halvade.hadoop.partitioners.*;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.util.Tool;
import be.ugent.intec.halvade.utils.Logger;
import be.ugent.intec.halvade.utils.HalvadeConf;
import be.ugent.intec.halvade.utils.Timer;
import org.seqdoop.hadoop_bam.BAMInputFormat;
import org.seqdoop.hadoop_bam.BAMOutputFormat;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import org.apache.hadoop.mapreduce.Job;

/**
 *
 * @author ddecap
 */
public class MapReduceRunner extends Configured implements Tool {
	protected final String DNA = " DNA job";
	protected HalvadeOptions halvadeOpts;

	@Override
	public int run(String[] strings) throws Exception {
		int ret = 0;
		try {
			Configuration halvadeConf = getConf();
			halvadeOpts = new HalvadeOptions();
			int optReturn = halvadeOpts.GetOptions(strings, halvadeConf);
			if (optReturn != 0)
				return optReturn;

			String halvadeDir = halvadeOpts.out + "bamfiles";
			ret = runHalvadeJob(halvadeConf, halvadeDir, HalvadeResourceManager.DNA);
			if (ret != 0) {
				Logger.DEBUG("Halvade align and filter job failed.");
				System.exit(-2);
			}
		} catch (IOException | ClassNotFoundException | IllegalArgumentException | IllegalStateException | InterruptedException | URISyntaxException e) {
			Logger.EXCEPTION(e);
		}
		return ret;
	}

	protected int runHalvadeJob(Configuration halvadeConf, String tmpOutDir, int jobType) throws IOException, URISyntaxException, InterruptedException, ClassNotFoundException {
		String pipeline = "DNA";
		HalvadeResourceManager.setJobResources(halvadeOpts, halvadeConf, jobType, false, halvadeOpts.useBamInput);

		HalvadeConf.setOutDir(halvadeConf, tmpOutDir);
		FileSystem outFs = FileSystem.get(new URI(tmpOutDir), halvadeConf);
		if (outFs.exists(new Path(tmpOutDir))) {
			Logger.INFO("The output directory \'" + tmpOutDir + "\' already exists.");
			Logger.INFO("ERROR: Please remove this directory before trying again.");
			System.exit(-2);
		}
		if (halvadeOpts.useBamInput)
			setHeaderFile(halvadeOpts.in, halvadeConf);
		Job halvadeJob = Job.getInstance(halvadeConf, "Halvade" + pipeline);
		halvadeJob.addCacheArchive(new URI(halvadeOpts.halvadeBinaries));
		halvadeJob.setJarByClass(be.ugent.intec.halvade.hadoop.mapreduce.HalvadeMapper.class);
		addInputFiles(halvadeOpts.in, halvadeConf, halvadeJob);
		FileOutputFormat.setOutputPath(halvadeJob, new Path(tmpOutDir));

		halvadeJob.setMapperClass(halvadeOpts.alignmentTools[halvadeOpts.aln]);
		halvadeJob.setReducerClass(be.ugent.intec.halvade.hadoop.mapreduce.filterReadsReducer.class);

		if (halvadeOpts.justAlign)
			halvadeJob.setNumReduceTasks(0);
		else
			halvadeJob.setNumReduceTasks(halvadeOpts.reduces);

		halvadeJob.setMapOutputKeyClass(ChromosomeRegion.class);
		halvadeJob.setMapOutputValueClass(SAMRecordWritable.class);
		halvadeJob.setInputFormatClass(HalvadeTextInputFormat.class);
		halvadeJob.setPartitionerClass(ChrRgPartitioner.class);
		halvadeJob.setSortComparatorClass(ChrRgSortComparator.class);
		halvadeJob.setGroupingComparatorClass(ChrRgGroupingComparator.class);
		halvadeJob.setOutputKeyClass(Text.class);
		halvadeJob.setOutputValueClass(BAMOutputFormat.class);

		if (halvadeOpts.useBamInput) {
			halvadeJob.setMapperClass(be.ugent.intec.halvade.hadoop.mapreduce.AlignedBamMapper.class);
			halvadeJob.setInputFormatClass(BAMInputFormat.class);
		}

		return runTimedJob(halvadeJob, "Halvade Job");
	}

	protected int runTimedJob(Job job, String jobname)
			throws IOException, InterruptedException, ClassNotFoundException {
		if (halvadeOpts.dryRun)
			return 0;
		Logger.DEBUG("Started " + jobname);
		Timer timer = new Timer();
		timer.start();
		int ret = job.waitForCompletion(true) ? 0 : 1;
		timer.stop();
		Logger.DEBUG("Finished " + jobname + " [runtime: " + timer.getFormattedElapsedTime() + "]");
		return ret;
	}

	protected void setHeaderFile(String input, Configuration conf) throws IOException, URISyntaxException {
		FileSystem fs = FileSystem.get(new URI(input), conf);
		String headerFile = null;
		if (fs.getFileStatus(new Path(input)).isDirectory()) {
			FileStatus[] files = fs.listStatus(new Path(input));
			if (files.length > 0)
				headerFile = files[0].getPath().toString();
		} else
			headerFile = input;
		if (headerFile != null)
			HalvadeConf.setHeaderFile(conf, headerFile);
	}

	protected void addInputFiles(String input, Configuration conf, Job job) throws URISyntaxException, IOException {
		FileSystem fs = FileSystem.get(new URI(input), conf);
		if (fs.getFileStatus(new Path(input)).isDirectory()) {
			// add every file in directory
			FileStatus[] files = fs.listStatus(new Path(input));
			for (FileStatus file : files) {
				if (!file.isDirectory()) {
					FileInputFormat.addInputPath(job, file.getPath());
				}
			}
		} else
			FileInputFormat.addInputPath(job, new Path(input));
	}

}
