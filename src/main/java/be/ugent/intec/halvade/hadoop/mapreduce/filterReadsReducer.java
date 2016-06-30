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
package be.ugent.intec.halvade.hadoop.mapreduce;

import be.ugent.intec.halvade.utils.SAMRecordIterator;
import org.seqdoop.hadoop_bam.SAMRecordWritable;
import be.ugent.intec.halvade.hadoop.datatypes.ChromosomeRegion;
import java.io.File;
import java.io.IOException;
import htsjdk.samtools.*;
import org.apache.hadoop.fs.FileSystem;
import be.ugent.intec.halvade.tools.PreprocessingTools;
import be.ugent.intec.halvade.tools.ProcessException;
import be.ugent.intec.halvade.tools.QualityException;
import be.ugent.intec.halvade.utils.ChromosomeRange;
import be.ugent.intec.halvade.utils.HalvadeFileUtils;
import be.ugent.intec.halvade.utils.HalvadeConf;
import be.ugent.intec.halvade.utils.Logger;
import java.net.URI;
import java.net.URISyntaxException;

/**
 *
 * @author ddecap
 */
public class filterReadsReducer extends HalvadeReducer {

    protected boolean isFirstAttempt;
    protected boolean filterDBsnp;
    protected boolean useUnifiedGenotyper;
    protected boolean redistribute;
    protected int newMaxQualScore = 60;
    protected int windows, cluster;
    protected double minFS, maxQD;
    protected int containers;
    protected int tasksLeft;
    protected String filterBedFile;

    @Override
    protected void reduce(ChromosomeRegion key, Iterable<SAMRecordWritable> values, Context context) throws IOException, InterruptedException {
        super.reduce(key, values, context);
        try {
            Logger.DEBUG("Processing key: " + key);
            // wrappers to call external programs
            PreprocessingTools tools = new PreprocessingTools(bin);
            tools.setContext(context);
            processAlignments(values, context, tools);
        } catch (URISyntaxException | QualityException | ProcessException ex) {
            Logger.EXCEPTION(ex);
            throw new InterruptedException(ex.getMessage());
        }
    }

    @Override
    protected void setup(Context context) throws IOException, InterruptedException {
        super.setup(context);
        isFirstAttempt = taskId.endsWith("_0");
        filterBedFile = HalvadeConf.getFilterBed(context.getConfiguration());
        filterDBsnp = HalvadeConf.getFilterDBSnp(context.getConfiguration());
        useUnifiedGenotyper = HalvadeConf.getUseUnifiedGenotyper(context.getConfiguration());
        redistribute = HalvadeConf.getRedistribute(context.getConfiguration());
        containers = HalvadeConf.getMapContainerCount(context.getConfiguration());
        tasksLeft = Integer.parseInt(context.getConfiguration().get("mapred.map.tasks")) - taskNr;
        // get task number: 
        if (redistribute && tasksLeft < containers) {
            threads = 6;
        }
    }

    protected void processAlignments(Iterable<SAMRecordWritable> values, Context context, PreprocessingTools tools)
            throws IOException, InterruptedException, URISyntaxException, QualityException {
        long startTime = System.currentTimeMillis();
        // temporary files
        String resultBam = tmpFileBase + ".bam";   
        boolean useElPrep = HalvadeConf.getUseElPrep(context.getConfiguration());
        ChromosomeRange r = new ChromosomeRange();
        SAMRecordIterator SAMit = new SAMRecordIterator(values.iterator(), header, r);        
        if(useElPrep && isFirstAttempt) 
            elPrepPreprocess(context, tools, SAMit, resultBam);
        else  {
            if(!isFirstAttempt) Logger.DEBUG("attempt " + taskId + ", preprocessing with Picard for smaller peak memory");
            PicardPreprocess(context, tools, SAMit, resultBam);
        }
        File result = new File(resultBam);
        if (result.exists()) {
        	bamFiles.add(resultBam);
        }
        long estimatedTime = System.currentTimeMillis() - startTime;
        Logger.DEBUG("total estimated time: " + estimatedTime / 1000);
    }

    protected void elPrepPreprocess(Context context, PreprocessingTools tools, SAMRecordIterator input, String output) throws InterruptedException, IOException, QualityException, URISyntaxException {
        String dictF = ref.substring(0, ref.lastIndexOf('.')) + ".dict";
        String rg = createReadGroupRecordString(RGID, RGLB, RGPL, RGPU, RGSM);
        String preSamOut = tmpFileBase + "-p1.sam";
        String samOut = tmpFileBase + "-p2.sam";
        String fCounts = tmpFileBase + "-features.count";

        outHeader = header.clone();
        outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        Logger.DEBUG("call elPrep");
        context.setStatus("call elPrep");
        int reads;
        if (keep) {
            reads = tools.callElPrep(preSamOut, samOut, inputIsBam ? null : rg, threads, input, outHeader, dictF);
        } else {
            reads = tools.streamElPrep(context, samOut, inputIsBam ? null : rg, threads, input, outHeader, dictF);
        }

        Logger.DEBUG(reads + " reads processed in elPrep");
        context.getCounter(HalvadeCounters.IN_PREP_READS).increment(reads);

        context.setStatus("convert SAM to BAM");
        Logger.DEBUG("convert SAM to BAM");
        tools.callSAMToBAM(samOut, output, threads);
        // remove temporary files
        HalvadeFileUtils.removeLocalFile(keep, preSamOut, context, HalvadeCounters.FOUT_GATK_TMP);
        HalvadeFileUtils.removeLocalFile(keep, samOut, context, HalvadeCounters.FOUT_GATK_TMP);
        HalvadeFileUtils.removeLocalFile(keep, fCounts);
    }

    protected void PicardPreprocess(Context context, PreprocessingTools tools, SAMRecordIterator input, String output) throws InterruptedException, QualityException, IOException, URISyntaxException {
        outHeader = header.clone();
        outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        // tmp files
        String tmpOut1 = tmpFileBase + "-p1.bam";
        String tmpOut2 = tmpFileBase + "-p2.bam";
        String tmpOut3 = tmpFileBase + "-p3.sam";
        String fCounts = tmpFileBase + "-features.count";
        String tmpMetrics = tmpFileBase + "-p3-metrics.txt";
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        if (!inputIsBam) {
            outHeader.addReadGroup(bamrg);
        }
        SAMFileWriter writer = factory.makeBAMWriter(outHeader, true, new File(tmpOut1));

        long startTime = System.currentTimeMillis();

        int count = 0;
        SAMRecord sam;
        while (input.hasNext()) {
            sam = input.next();
            writer.addAlignment(sam);
            count++;
        }
        int reads = input.getCount();
        writer.close();

        context.getCounter(HalvadeCounters.IN_PREP_READS).increment(reads);
        long estimatedTime = System.currentTimeMillis() - startTime;
        context.getCounter(HalvadeCounters.TIME_HADOOP_SAMTOBAM).increment(estimatedTime);
        Logger.DEBUG("time writing " + count + " records to disk: " + estimatedTime / 1000);

        Logger.DEBUG("clean sam");
        context.setStatus("clean sam");
        tools.runCleanSam(tmpOut1, tmpOut2);
        Logger.DEBUG("mark duplicates");
        context.setStatus("mark duplicates");
        tools.runMarkDuplicates(tmpOut2, tmpOut3, tmpMetrics);

        if (!inputIsBam) {
            Logger.DEBUG("add read-group");
            context.setStatus("add read-group");
            tools.runAddOrReplaceReadGroups(tmpOut3, output, RGID, RGLB, RGPL, RGPU, RGSM);
        } else {
            context.setStatus("convert SAM to BAM");
            Logger.DEBUG("convert SAM to BAM");
            tools.callSAMToBAM(tmpOut3, output, threads);
        }

        estimatedTime = System.currentTimeMillis() - startTime;
        Logger.DEBUG("estimated time: " + estimatedTime / 1000);

        // remove all temporary files now!
        HalvadeFileUtils.removeLocalFile(keep, tmpMetrics, context, HalvadeCounters.FOUT_GATK_TMP);
        HalvadeFileUtils.removeLocalFile(keep, tmpOut1, context, HalvadeCounters.FOUT_GATK_TMP);
        HalvadeFileUtils.removeLocalFile(keep, tmpOut2, context, HalvadeCounters.FOUT_GATK_TMP);
        HalvadeFileUtils.removeLocalFile(keep, tmpOut3, context, HalvadeCounters.FOUT_GATK_TMP);
        HalvadeFileUtils.removeLocalFile(keep, fCounts);
    }

    protected String makeRegionFile(Context context, ChromosomeRange r, PreprocessingTools tools, String region) throws URISyntaxException, IOException, InterruptedException {
        // if exome dont do but for exome filter on exomeBedFile
        if (filterBedFile == null) {
            r.writeToPicardRegionFile(region);
        } else {
            String exomebed = tmpFileBase + "exome.bed";
            if (filterBedFile.endsWith(".gz")) {
                exomebed += ".gz";
            }
            HalvadeFileUtils.downloadFileFromHDFS(context, FileSystem.get(new URI(filterBedFile), context.getConfiguration()),
                    filterBedFile, exomebed);
            if (exomebed.endsWith(".gz")) {
                exomebed = HalvadeFileUtils.Unzip(exomebed);
            }
            region = tools.filterExomeBed(exomebed, r);
            if (region == null) {
                Logger.DEBUG("empty region file, no vcf results!!");
                return null;
            }
            HalvadeFileUtils.removeLocalFile(keep, exomebed);
        }
        return region;
    }

}
