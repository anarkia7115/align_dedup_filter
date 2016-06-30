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


import be.ugent.intec.halvade.utils.ChromosomeSplitter;
import be.ugent.intec.halvade.utils.Logger;
import be.ugent.intec.halvade.utils.HalvadeConf;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.text.DecimalFormat;
import java.util.Enumeration;
import java.util.Properties;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FSDataInputStream;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.mapreduce.Mapper;

/**
 *
 * @author ddecap
 */
public class HalvadeOptions {

    public Options options = new Options();
    public String in;
    public String out;
    public String ref;
    public String STARGenome = null;
	public String reduceTool = null;
    public String java = null;
    public String gff = null;
    public String tmpDir = "/tmp/halvade/";
    public String localRefDir = null;
    public String sites;
    public int nodes, vcores;
    public double mem;
    public int maps = 1, reduces = 1, mthreads = 1, rthreads = 1;
    public String[] hdfsSites;
    public boolean paired = true;
    public int aln = 0;
    public Class<? extends Mapper>[] alignmentTools = new Class[]{
        be.ugent.intec.halvade.hadoop.mapreduce.BWAAlnMapper.class,
        be.ugent.intec.halvade.hadoop.mapreduce.BWAMemMapper.class,
        be.ugent.intec.halvade.hadoop.mapreduce.Bowtie2Mapper.class,
        be.ugent.intec.halvade.hadoop.mapreduce.Cushaw2Mapper.class
    };
    public boolean justCombine = false;
    public boolean filterDBSnp = false;
    public boolean useGenotyper = true;
    public String RGID = "GROUP1";
    public String RGLB = "LIB1";
    public String RGPL = "ILLUMINA";
    public String RGPU = "UNIT1";
    public String RGSM = "SAMPLE1";
    public boolean useElPrep = true;
    public boolean keepFiles = false;
    public int stand_call_conf = -1;
    public int stand_emit_conf = -1;
    public SAMSequenceDictionary dict;
    public String chr = null;
    public int reducerContainersPerNode = -1;
    public int mapContainersPerNode = -1;
    public boolean justAlign = false;
    public String bedFile = null;
    public String filterBed = null;
    public String bedRegion = null;
    public double coverage = -1.0;
    public String halvadeBinaries;
    public String bin;
    public boolean hdfspath = true;
    public boolean combineVcf = false;
    public boolean dryRun = false;
    public boolean keepChrSplitPairs = true;
    public boolean startupJob = true;
    public boolean rnaPipeline = false;
    public boolean reportAll = false;
    public boolean useBamInput = false;
    public boolean setMapContainers = true, setReduceContainers = true;
    public boolean redistribute = false;
    public boolean smtEnabled = false;
    public int overrideMem = -1;
    
    protected DecimalFormat onedec;
    protected static final double REDUCE_TASKS_FACTOR = 1.68 * 15;
    protected static final double DEFAULT_COVERAGE = 50;
    protected static final double DEFAULT_COVERAGE_SIZE = 86;
    protected static final String DICT_SUFFIX = ".dict";

    public int GetOptions(String[] args, Configuration hConf) throws IOException, URISyntaxException {
        try {
            boolean result = parseArguments(args, hConf);
            if (!result) {
                HelpFormatter formatter = new HelpFormatter();
                formatter.setWidth(80);
                formatter.printHelp("hadoop jar HalvadeWithLibs.jar -I <IN> -O <OUT> "
                        + "-R <REF> -D <SITES> -B <BIN> -nodes <nodes> -mem <mem> -vcores <cores> [options]", options);
                return 1;
            }
            onedec = new DecimalFormat("###0.0");
            // add parameters to configuration:
            if (localRefDir == null) {
                localRefDir = tmpDir;
            }
            HalvadeConf.setScratchTempDir(hConf, tmpDir);
            HalvadeConf.setRefDirOnScratch(hConf, localRefDir);
            HalvadeConf.setRefOnHDFS(hConf, ref);
            if (STARGenome != null) {
                HalvadeConf.setStarDirOnHDFS(hConf, STARGenome);
            }
            HalvadeConf.setKnownSitesOnHDFS(hConf, hdfsSites);
            HalvadeConf.setIsPaired(hConf, paired);
            HalvadeConf.setIsRNA(hConf, rnaPipeline);
            if (bedFile != null) {            
                HalvadeConf.setBed(hConf, bedFile);
            }
            if (filterBed != null) {            
                HalvadeConf.setFilterBed(hConf, filterBed);
            }
            HalvadeConf.setInputIsBam(hConf, useBamInput);
            HalvadeConf.setOutDir(hConf, out);
            HalvadeConf.setKeepFiles(hConf, keepFiles);
            HalvadeConf.setFilterDBSnp(hConf, filterDBSnp);
            HalvadeConf.clearTaskFiles(hConf);
            HalvadeConf.setUseElPrep(hConf, useElPrep);
            HalvadeConf.setUseUnifiedGenotyper(hConf, useGenotyper);
            HalvadeConf.setRedistribute(hConf, redistribute);
            HalvadeConf.setReadGroup(hConf, "ID:" + RGID + " LB:" + RGLB + " PL:" + RGPL + " PU:" + RGPU + " SM:" + RGSM);
            HalvadeConf.setkeepChrSplitPairs(hConf, keepChrSplitPairs);
            if (STARGenome != null) {
                HalvadeConf.setStarDirPass2HDFS(hConf, out);
            }

            if (chr != null) {
                HalvadeConf.setChrList(hConf, chr);
            }
            if (java != null) {
                HalvadeConf.setJava(hConf, java);
            }
            if (gff != null) {
                HalvadeConf.setGff(hConf, gff);
            }

            if (stand_call_conf > 0) {
                HalvadeConf.setSCC(hConf, stand_call_conf);
            }
            if (stand_emit_conf > 0) {
                HalvadeConf.setSEC(hConf, stand_emit_conf);
            }

            parseDictFile(hConf);
            double inputSize = getInputSize(in, hConf);
            if (coverage == -1.0) {
                coverage = Math.max(1.0, DEFAULT_COVERAGE * (inputSize / DEFAULT_COVERAGE_SIZE));
            }
            Logger.DEBUG("Estimated coverage: " + roundOneDecimal(coverage));
            // set a minimum first where the real amount is based on
            reduces = (int) (coverage * REDUCE_TASKS_FACTOR);
            Logger.DEBUG("estimated # reducers: " + reduces);
            ChromosomeSplitter splitter;
            if(bedFile != null && hdfspath) {
                // Check the local path exist ...
                //File file = new File(tmpDir);
                //if (!file.exists()  && ! file.isDirectory()) 
                //	Logger.DEBUG("The tmpDir does not exist before download the bedfile.");
                String newBedFile = tmpDir + bedFile.substring(bedFile.lastIndexOf("/")+1);
                FileSystem fs = FileSystem.get(new URI(bedFile), hConf);
                int beddownload = downloadFileFromHDFS(fs, bedFile, newBedFile);
                if (beddownload == 0) {
                	Logger.DEBUG("Download the hdfs bedfile to local disk successfully.");
                } else {
                	Logger.DEBUG("Error occurs in downloading the hdfs bedfile to local disk !");
                }
                splitter = new ChromosomeSplitter(dict, newBedFile, reduces);
            } else if (bedFile != null) {
            	splitter = new ChromosomeSplitter(dict, bedFile, reduces);
            } else {
                splitter = new ChromosomeSplitter(dict, reduces);
            }
            String bedRegions = out + "HalvadeRegions.bed";
            splitter.exportSplitter(bedRegions, hConf);
            reduces = splitter.getRegionCount();
            Logger.DEBUG("actual # reducers: " + reduces);
            HalvadeConf.setBedRegions(hConf, bedRegions);

        } catch (ParseException e) {
            Logger.DEBUG(e.getMessage());
            HelpFormatter formatter = new HelpFormatter();
            formatter.setWidth(80);
            formatter.printHelp("hadoop jar HalvadeWithLibs.jar -I <input> -O <output> "
                    + "-R <ref> -D <dbsnp> -B <bin> -nodes <nodes> -mem <mem> -vcores <cores> [options]", options);
            return 1;
        }
        return 0;
    }

    protected String roundOneDecimal(double val) {
        return onedec.format(val);
    }

    protected double getInputSize(String input, Configuration conf) throws URISyntaxException, IOException {
        double size = 0;
        FileSystem fs = FileSystem.get(new URI(input), conf);
        if (fs.getFileStatus(new Path(input)).isDirectory()) {
            // add every file in directory
            FileStatus[] files = fs.listStatus(new Path(input));
            for (FileStatus file : files) {
                if (!file.isDirectory()) {
                    size += file.getLen();
                }
            }
        } else {
            size += fs.getFileStatus(new Path(input)).getLen();
        }
        return (size / (1024 * 1024 * 1024));
    }


    protected void parseDictFile(Configuration conf) {
        be.ugent.intec.halvade.utils.Logger.DEBUG("parsing dictionary " + ref + DICT_SUFFIX);
        try {
            FileSystem fs = FileSystem.get(new URI(ref + DICT_SUFFIX), conf);
            FSDataInputStream stream = fs.open(new Path(ref + DICT_SUFFIX));
            String line = getLine(stream); // header
            dict = new SAMSequenceDictionary();
            line = getLine(stream);
            while (line != null) {
                String[] lineData = line.split("\\s+");
                String seqName = lineData[1].substring(lineData[1].indexOf(':') + 1);
                int seqLength = 0;
                try {
                    seqLength = Integer.parseInt(lineData[2].substring(lineData[2].indexOf(':') + 1));
                } catch (NumberFormatException ex) {
                    be.ugent.intec.halvade.utils.Logger.EXCEPTION(ex);
                }
                SAMSequenceRecord seq = new SAMSequenceRecord(seqName, seqLength);
//                Logger.DEBUG("name: " + seq.getSequenceName() + " length: " + seq.getSequenceLength());
                dict.addSequence(seq);
                line = getLine(stream);
            }
            HalvadeConf.setSequenceDictionary(conf, dict);
        } catch (URISyntaxException | IOException ex) {
            be.ugent.intec.halvade.utils.Logger.EXCEPTION(ex);
        }

    }

    protected String getLine(FSDataInputStream stream) throws IOException {
        String tmp = "";
        try {
            char c = (char) stream.readByte();
            while (c != '\n') {
                tmp = tmp + c;
                c = (char) stream.readByte();
            }
            return tmp;
        } catch (EOFException ex) {
            // reached end of file, return null;
            return null;
        }
    }
    
    /**
     * @return returns 0 if successfull, -1 if filesize is incorrect and -2 if an exception occurred
     */
    private static int privateDownloadFileFromHDFS(FileSystem fs, String from, String to) {
        try {
            // check if file is present on local scratch
            File f = new File(to);
            if(!f.exists()) {
                Logger.DEBUG("attempting download of \"" + to + "\"");
                fs.copyToLocalFile(new Path(from), new Path(to));
            } else {
                // check if filesize is correct
                if(fs.getFileStatus(new Path(from)).getLen() != f.length()) {
                    // incorrect filesize, remove and download again
                    Logger.DEBUG("incorrect filesize: " + f.length() + " =/= " + 
                            fs.getFileStatus(new Path(from)).getLen());
                    f.delete();
                    fs.copyToLocalFile(new Path(from), new Path(to));
            
                } else {
                    Logger.DEBUG("file \"" + to + "\" exists");
                }
            }
            if(fs.getFileStatus(new Path(from)).getLen() != f.length())
                return -1;
            else
                return 0;
        } catch (IOException ex) {
            Logger.DEBUG("failed to download " + from + " from HDFS: " + ex.getLocalizedMessage());
            Logger.EXCEPTION(ex);
            return -2;
        }
    }
    private static int attemptDownloadFileFromHDFS(FileSystem fs, String from, String to, int tries) throws IOException {
        if(from.equalsIgnoreCase(to)) return 0;
        int val = privateDownloadFileFromHDFS(fs, from, to);
        int try_ = 1;
        while (val != 0 && try_ < tries) {
            val = privateDownloadFileFromHDFS(fs, from, to);
            try_++;
        }
        if(val == 0)
            Logger.DEBUG(from + " downloaded");
        else {
            Logger.DEBUG(from + " failed to download");   
            throw new IOException();
        }
        return val;            
    }
    public static int downloadFileFromHDFS(FileSystem fs, String from, String to) throws IOException {
        return attemptDownloadFileFromHDFS(fs, from, to, 3);            
    }
    
    protected void createOptions() {
        OptionBuilder.withArgName("input");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .isRequired(true);
		OptionBuilder
                .withDescription("Input directory on hdfs containing fastq files.");
		Option optIn = OptionBuilder
                .create("I");
        OptionBuilder.withArgName("output");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .isRequired(true);
		OptionBuilder
                .withDescription("Output directory on hdfs.");
		Option optOut = OptionBuilder
                .create("O");
        OptionBuilder.withArgName("bin.tar.gz");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .isRequired(true);
		OptionBuilder
                .withDescription("The tarred file containing all binary files located on HDFS.");
		Option optBin = OptionBuilder
                .create("B");
        OptionBuilder.withArgName("reference");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .isRequired(true);
		OptionBuilder
                .withDescription("Name of the fasta file name of the reference (without extension) on HDFS. Make sure the BWA index has the same prefix.");
		Option optRef = OptionBuilder
                .create("R");
        OptionBuilder.withArgName("nodes");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .isRequired(true);
		OptionBuilder
                .withDescription("Sets the number of nodes in this cluster.");
		Option optNodes = OptionBuilder
                .create("nodes");
        OptionBuilder.withArgName("cores");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .isRequired(true);
		OptionBuilder
                .withDescription("Sets the available cpu cores per node in this cluster.");
		Option optVcores = OptionBuilder
                .create("vcores");
        OptionBuilder.withArgName("gb");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .isRequired(true);
		OptionBuilder
                .withDescription("Sets the available memory [in GB] per node in this cluster.");
		Option optMem = OptionBuilder
                .create("mem");
        OptionBuilder.withArgName("gb");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Overrides the maximum container memory [in GB].");
		Option optRmem = OptionBuilder
                .create("refmem");
        OptionBuilder.withArgName("snpDBa,snpDBb");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .isRequired(true);
		OptionBuilder
                .withDescription("Name of snpDB files for the genome on HDFS. If multiple separate with \',\'.");
		Option optSites = OptionBuilder
                .create("D");
        OptionBuilder.withArgName("stargenome");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Directory on HDFS containing all STAR genome files");
		Option optStarGenome = OptionBuilder
                .create("SG");
        OptionBuilder.withArgName("dir");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets the location for temporary files on every node [/tmp/halvade/].");
		Option optTmp = OptionBuilder
                .create("tmp");
        OptionBuilder.withArgName("dir");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets the folder containing all the reference files for BWA or STAR and GATK on every node [tmp directory].");
		Option optrefdir = OptionBuilder
                .create("refdir");
        OptionBuilder.withArgName("java");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Set location of java binary to use [must be 1.7+].");
		Option optJava = OptionBuilder
                .create("Java");
        OptionBuilder.withArgName("Set the reducer tool for snp/indel analysis.");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("We support the snp/indel analysis for GATK and bcftools currently.");
		Option optRT = OptionBuilder
                .create("reduceTool");
        OptionBuilder.withArgName("RGID");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("sets the RGID for the read-group.");
		// "ID:" + RGID + " LB:" + RGLB + " PL:" + RGPL + " PU:" + RGPU + " SM:" + RGSM
        Option optID = OptionBuilder
                .create("id");
        OptionBuilder.withArgName("RGLB");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("sets the RGLB for the read-group.");
		Option optLB = OptionBuilder
                .create("lb");
        OptionBuilder.withArgName("RGPL");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("sets the RGPL for the read-group.");
		Option optPL = OptionBuilder
                .create("pl");
        OptionBuilder.withArgName("RGPU");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("sets the RGPU for the read-group.");
		Option optPU = OptionBuilder
                .create("pu");
        OptionBuilder.withArgName("RGSM");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("sets the RGSM for the read-group.");
		Option optSM = OptionBuilder
                .create("sm");
        OptionBuilder.withArgName("coverage");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets the coverage to better distribute the tasks.");
		Option optCov = OptionBuilder
                .create("cov");
        OptionBuilder.withArgName("scc");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets stand_call_conf for gatk Variant Caller.");
		Option optScc = OptionBuilder
                .create("scc");
        OptionBuilder.withArgName("sec");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets stand_emit_conf for gatk Variant Caller.");
		Option optSec = OptionBuilder
                .create("sec");
        OptionBuilder.withArgName("bedfile");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets the bed file containing relevant (Genes) regions which "
                        + "will be used to split the genome into genomic regions.");
		Option optBed = OptionBuilder
                .create("bed");
        OptionBuilder.withArgName("bedfile");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets the bed file containing relevant (Exome) regions which "
                        + " will be used to filter in the GATK steps.");
		Option optFilterBed = OptionBuilder
                .create("fbed");
        OptionBuilder.withArgName("gff");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets the gff file to be used with HTSeq-Count. This is required to run HTSeq-Count.");
		Option optGff = OptionBuilder
                .create("gff");
        OptionBuilder.withArgName("tasks");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Overrides the number of map tasks running simultaneously on each node. ");
		Option optMpn = OptionBuilder
                .create("mpn");
        OptionBuilder.withArgName("tasks");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Overrides the number of reduce tasks running simultaneously on each node. ");
		Option optRpn = OptionBuilder
                .create("rpn");
        OptionBuilder.withLongOpt("custom_args");
		OptionBuilder
                .withArgName("tool=args");
		OptionBuilder
                .hasArgs(2);
		OptionBuilder
                .withValueSeparator();
		OptionBuilder
                .withDescription("Adds custom arguments for a tool. If a module in a tool is used, add the name after an underscore. "
                        + "Possible values: " + getProgramNames());
		Option optCustomArgs = OptionBuilder
                .create("CA");
        OptionBuilder.withArgName("num");
		OptionBuilder
                .hasArg();
		OptionBuilder
                .withDescription("Sets the aligner used in Halvade. Possible values are 0 (bwa aln+sampe)[default], 1 (bwa mem), 2 (bowtie2), 3 (cushaw2).");
		Option optAln = OptionBuilder
                .create("aln");

        OptionBuilder.withDescription("Sets the input files to single reads [default is paired-end reads].");
		//flags
        Option optSingle = OptionBuilder
                .create("single");
        OptionBuilder.withDescription("Just Combines the vcf on HDFS [out dir] and doesn't run the hadoop job.");
		Option optCombine = OptionBuilder
                .create("combine");
        OptionBuilder.withDescription("Uses Picard to preprocess the data for GATK.");
		Option optPp = OptionBuilder
                .create("Picard");
        OptionBuilder.withDescription("Use Bedtools to select the needed interval of dbsnp.");
		Option optFilterDBsnp = OptionBuilder
                .create("filter_dbsnp");
        OptionBuilder.withDescription("Only align the reads.");
		Option optJustAlign = OptionBuilder
                .create("justalign");
        OptionBuilder.withDescription("Enable simultaneous multithreading.");
		Option optSmt = OptionBuilder
                .create("smt");
        OptionBuilder.withDescription("Keep intermediate files.");
		Option optKeep = OptionBuilder
                .create("keep");
        OptionBuilder.withDescription("Use HaplotypeCaller instead of UnifiedGenotyper for Variant Detection.");
		Option optHap = OptionBuilder
                .create("hc");
        OptionBuilder.withDescription("Run the RNA Best Practices pipeline by Broad [default is DNA pipeline]. SG needs to be set for this.");
		Option optRna = OptionBuilder
                .create("rna");
        OptionBuilder.withDescription("Execute a dryrun, will calculate task size, split for regions etc, but not execute the MapReduce job.");
		Option optDry = OptionBuilder
                .create("dryrun");
        OptionBuilder.withDescription("Drop all paired-end reads where the pairs are aligned to different chromosomes.");
		Option optDrop = OptionBuilder
                .create("drop");
        OptionBuilder.withDescription("Reports all variants at the same location when combining variants.");
		Option optReportAll = OptionBuilder
                .create("report_all");
        OptionBuilder.withDescription("Uses aligned bam as input files instead of unaligned fastq files.");
		Option optBamIn = OptionBuilder
                .create("bam");
        OptionBuilder.withDescription("This will enable Halvade to redistribute resources when possible when not all containers are used.");
		Option optRedis = OptionBuilder
                .create("redistribute");

        options.addOption(optIn);
        options.addOption(optOut);
        options.addOption(optRef);
        options.addOption(optSites);
        options.addOption(optBin);
        options.addOption(optTmp);
        options.addOption(optrefdir);
        options.addOption(optSingle);
        options.addOption(optAln);
        options.addOption(optRT);
        options.addOption(optID);
        options.addOption(optLB);
        options.addOption(optPL);
        options.addOption(optPU);
        options.addOption(optSM);
        options.addOption(optPp);
        options.addOption(optFilterDBsnp);
        options.addOption(optHap);
        options.addOption(optScc);
        options.addOption(optSec);
        options.addOption(optBed);
        options.addOption(optFilterBed);
        options.addOption(optJava);
        options.addOption(optCombine);
        options.addOption(optNodes);
        options.addOption(optVcores);
        options.addOption(optMem);
        options.addOption(optKeep);
        options.addOption(optJustAlign);
        options.addOption(optCov);
        options.addOption(optMpn);
        options.addOption(optGff);
        options.addOption(optRpn);
        options.addOption(optDry);
        options.addOption(optDrop);
        options.addOption(optReportAll);
        options.addOption(optSmt);
        options.addOption(optRna);
        options.addOption(optStarGenome);
        options.addOption(optBamIn);
        options.addOption(optCustomArgs);
        options.addOption(optRedis);
        options.addOption(optRmem);
    }

    protected boolean parseArguments(String[] args, Configuration halvadeConf) throws ParseException {
        createOptions();
        CommandLineParser parser = new GnuParser();
        CommandLine line = parser.parse(options, args);

        in = line.getOptionValue("I");
        out = line.getOptionValue("O");
        if(!out.endsWith("/")) out += "/";
        ref = line.getOptionValue("R");
        sites = line.getOptionValue("D");
        halvadeBinaries = line.getOptionValue("B");
        hdfsSites = sites.split(",");
        if (line.hasOption("bam")) {
            useBamInput = true;
        }
        if (line.hasOption("rna")) {
            rnaPipeline = true;
            if (line.hasOption("SG"))
                STARGenome = line.getOptionValue("SG");
            if(!useBamInput && STARGenome == null) {
                throw new ParseException("the '-rna' option requires -SG, pointing to the location of the STAR reference directory if alignment with STAR aligner is required (fastq).");
            }
        }

        if (line.hasOption("tmp")) {
            tmpDir = line.getOptionValue("tmp");
        }
        if (line.hasOption("refdir")) {
            localRefDir = line.getOptionValue("refdir");
        }
        if (line.hasOption("RT")) {
        	reduceTool = line.getOptionValue("RT");
        }
        if (line.hasOption("nodes")) {
            nodes = Integer.parseInt(line.getOptionValue("nodes"));
        }
        if (line.hasOption("vcores")) {
            vcores = Integer.parseInt(line.getOptionValue("vcores"));
        }
        if (line.hasOption("smt")) {
            smtEnabled = true;
        }
        if (line.hasOption("mem")) {
            mem = Double.parseDouble(line.getOptionValue("mem"));
        }
        if (line.hasOption("mpn")) {
            setMapContainers = false;
            mapContainersPerNode = Integer.parseInt(line.getOptionValue("mpn"));
        }
        if (line.hasOption("rpn")) {
            setReduceContainers = false;
            reducerContainersPerNode = Integer.parseInt(line.getOptionValue("rpn"));
        }
        if (line.hasOption("refmem")) {
            overrideMem = Integer.parseInt(line.getOptionValue("refmem")) * 1024;
        }
        if (line.hasOption("scc")) {
            stand_call_conf = Integer.parseInt(line.getOptionValue("scc"));
        }
        if (line.hasOption("sec")) {
            stand_emit_conf = Integer.parseInt(line.getOptionValue("sec"));
        }
        if (line.hasOption("report_all")) {
            reportAll = true;
        }
        if (line.hasOption("keep")) {
            keepFiles = true;
        }
        if (line.hasOption("redistribute")) {
            redistribute = true;
        }
        if (line.hasOption("single")) {
            paired = false;
        }
        if (line.hasOption("justalign")) {
            justAlign = true;
            combineVcf = false;
        }
        if (line.hasOption("aln")) {
            aln = Integer.parseInt(line.getOptionValue("aln"));
            if(aln < 0 || aln > 3)
                aln = 0; // default value
        }
        if (line.hasOption("Java")) {
            java = line.getOptionValue("Java");
        }
        if (line.hasOption("gff")) {
            gff = line.getOptionValue("gff");
        }
        if (line.hasOption("dryrun")) {
            dryRun = true;
            combineVcf = false;
        }
        if (line.hasOption("drop")) {
            keepChrSplitPairs = false;
        }
        if (line.hasOption("cov")) {
            coverage = Integer.parseInt(line.getOptionValue("cov"));
        }
        if (line.hasOption("combine")) {
            justCombine = true;
            combineVcf = true;
        }
        if (line.hasOption("filter_dbsnp")) {
            filterDBSnp = true;
        }
        if (line.hasOption("hc")) {
            useGenotyper = false;
        }
        if (line.hasOption("Picard")) {
            useElPrep = false;
        }
        if (line.hasOption("id")) {
            RGID = line.getOptionValue("id");
        }
        if (line.hasOption("lb")) {
            RGLB = line.getOptionValue("lb");
        }
        if (line.hasOption("pl")) {
            RGPL = line.getOptionValue("pl");
        }
        if (line.hasOption("pu")) {
            RGPU = line.getOptionValue("pu");
        }
        if (line.hasOption("sm")) {
            RGSM = line.getOptionValue("sm");
        }
        if (line.hasOption("bed")) {
            bedFile = line.getOptionValue("bed");
        }
        if (line.hasOption("fbed")) {
            filterBed = line.getOptionValue("fbed");
        }

        if (line.hasOption("CA")) {
            Properties props = line.getOptionProperties("CA");
            Enumeration names = props.propertyNames();
            while (names.hasMoreElements()) {
                String name = (String) names.nextElement();
                addCustomArguments(halvadeConf, name, props.getProperty(name));
            }
        }
        return true;
    }

    protected String[] programNames = {
        "bwa_aln", "bwa_mem", "bwa_sampe",
        "star",
        "elprep",
        "samtools_view",
        "bedtools_bdsnp", "bedtools_exome",
        "picard_buildbamindex", "picard_addorreplacereadgroup", "picard_markduplicates", "picard_cleansam",
        "gatk_realignertargetcreator", "gatk_indelrealigner", "gatk_baserecalibrator", "gatk_printreads", "gatk_combinevariants",
        "gatk_variantcaller", "gatk_variantannotator", "gatk_variantfiltration", "gatk_splitncigarreads"};

    protected String getProgramNames() {
        String names = programNames[0];
        for (int i = 1; i < programNames.length; i++) {
            names += ", " + programNames[i];
        }
        return names;
    }

    protected void addCustomArguments(Configuration halvadeConf, String name, String property) throws ParseException {
        boolean found = false;
        int i = 0;
        while (!found && i < programNames.length) {
            if (name.equalsIgnoreCase(programNames[i])) {
                found = true;
            }
            i++;
        }
        if (found) {
            String[] split = name.split("_");
            String program = split[0];
            String tool = "";
            if (split.length > 1) {
                tool = split[1];
            }
            HalvadeConf.setCustomArgs(halvadeConf, program, tool, property);
            Logger.DEBUG("Custom arguments for " + name + ": \"" + property + "\"");
        } else {
            Logger.DEBUG("Unknown program: " + name + ", custom arguments [" + property + "] ignored");
            throw new ParseException("Program " + name + " not found, please use a valid name.");
        }
    }
}
