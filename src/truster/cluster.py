#! /bin/python3
import subprocess
import os
from .jobHandler import *
import pysam
from .bcolors import Bcolors

class Cluster:
    def __init__(self, cluster_name, tsv, logfile):
        self.cluster_name = cluster_name
        self.tsv = tsv
        self.outdirs = {"tsv_to_bam" : None, "filter_UMIs" : None, "bam_to_fastq" : None, "concatenate_lanes" : None, "map_cluster" : None, "TE_count_unique" : None, "TE_count" : None}
        self.logfile = logfile

    def tsv_to_bam(self, sample_id, bam, outdir, slurm=None, modules=None, remote_bam_flag = False, remote_bam_path = None, remote_bam_tmpdir = None, dry_run = False):
        if not os.path.exists("tsv_to_bam_scripts"):
            os.makedirs("tsv_to_bam_scripts", exist_ok=True)
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        with open(self.logfile, "a") as log:
            try:
                if not remote_bam_flag:
                    cmd = ["subset-bam", "--bam", bam, "--cell-barcodes", self.tsv, "--out-bam", os.path.join(outdir, (self.cluster_name + ".bam"))]
                    
                    result = run_instruction(cmd = cmd, fun = "tsv_to_bam", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "tsv_to_bam", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                    exit_code = result[1]
                else: # We have a remote bam file
                    # Let's copy it to its tmp dir using iget (bam path was already fixed in experiment.py tsv_to_bam_clusters()) 
                    cmd = ["iget", remote_bam_path, remote_bam_tmpdir, " || exit 2; ", "subset-bam", "--bam", bam, "--cell-barcodes", self.tsv, "--out-bam", os.path.join(outdir, (self.cluster_name + ".bam"))]
                    result = run_instruction(cmd = cmd, fun = "tsv_to_bam", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "tsv_to_bam", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                    exit_code = result[1]
                if exit_code == 0:
                    self.outdirs["tsv_to_bam"] = outdir
                return result
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def filter_UMIs(self, sample_id, inbam, outdir, slurm=None, modules=None, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("filter_UMIs_scripts"):
                    os.makedirs("filter_UMIs_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                outbam = os.path.join(outdir, (self.cluster_name + "_filtered.bam"))
                
                cwd = os.path.dirname(os.path.realpath(__file__))
                
                cmd = ["samtools view","-b", "-o", outbam, "-F", "0x400", inbam] # Remove PCR duplicates as annotated by cellranger
                result = run_instruction(cmd = cmd, fun = "filter_UMIs", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "filter_UMIs", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs["filter_UMIs"] = outdir
                return result

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def bam_to_fastq(self, sample_id, bam, outdir, slurm=None, modules=None, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("bam_to_fastq_scripts"):
                    os.makedirs("bam_to_fastq_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                cmd = ["bamtofastq", bam, (outdir + "/" + self.cluster_name)]
                result = run_instruction(cmd = cmd, fun = "bam_to_fastq", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "bam_to_fastq", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs["bam_to_fastq"] = outdir
                return result
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def concatenate_lanes(self, sample_id, indir, outdir, slurm=None, modules=None, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("concatenate_lanes_scripts"):
                    os.makedirs("concatenate_lanes_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)

                # Initialize the list of files to concatenate
                files_to_concatenate = []
                library_names = []
                # Walk through the files in the indir
                for root, subdirs, files in os.walk(os.path.join(indir, self.cluster_name)):
                    # If there is more than one subdirectory, we could be dealing with different library types
                    if len(subdirs) > 1: 
                        log.write("WARNING: I'm assuming library type 0 is gene expression. Please check I'm right. You have more than one library folder from sample " + sample_id + ", in: " + indir)
                    # The subdirectories have the form of sampleid_0_1_ID
                    for subdir in subdirs:
                        # We use sample_id to split as sample_id might contain underscores
                        library_type = subdir.split("_")[-3:] # Library ID, GEM group, Flowcell ID
                        library_type = library_type[:-1] # We dont care about the flowcell ID
                        if library_type == ["0", "1"]: # Gene expression
                            # For each of the gene expression libraries found in this sample, we walk through the files
                            for root_in_subdir, subdirs_in_subdir, files_in_subdir in os.walk(os.path.join(indir, self.cluster_name, subdir)):
                                for file in files_in_subdir:
                                    if(file.endswith("R2_001.fastq.gz")): # If it's a sequence file, we want to concatenate
                                        files_to_concatenate.append(os.path.join(indir, self.cluster_name, root_in_subdir, file))
                                        library_names.append(subdir)
                
                cwd = os.path.dirname(os.path.realpath(__file__))
                
                fastq_out = os.path.join(outdir, (self.cluster_name + "_R2.fastq.gz"))
                # We don't want to keep appending to an existing file...
                if os.path.exists(fastq_out):
                    print(f"Output file {fastq_out} exists. Please delete and try again.")
                    msg = f"Output file for concatenate_lanes {fastq_out} exists. Please delete and try again."
                    log.write(msg)
                    return("", 2) # Return error
                    
                cmd = ["cat", " ".join(files_to_concatenate), ">", fastq_out]
                if slurm is not None: # if we can use slurm, we do
                    result = run_instruction(cmd = cmd, fun = "concatenate_lanes", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "concatenate_lanes", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                    exit_code = result[1]
                    if exit_code == 0:
                        self.outdirs["concatenate_lanes"] = outdir
                else: # otherwise we need to reformulate the cat for subprocess.call
                    cmd = [fastq.replace(" ", "\ ") for fastq in files_to_concatenate]
                    cmd.insert(0, "cat")
                    with open(fastq_out, "w") as fout:
                        result = subprocess.call(cmd, shell = False, stdout=fout, universal_newlines=True)        
                    # If something went wrong, we return False
                    if result == 0:
                        self.outdirs["concatenate_lanes"] = outdir
                return result
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def map_cluster(self, sample_id, fastq_dir, outdir, gene_gtf, star_index, RAM, out_tmp_dir = None, unique=False, slurm=None, modules=None, snic_tmp = False, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("map_cluster_scripts"):
                    os.makedirs("map_cluster_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
    
                cmd = ["STAR", "--runThreadN", str(slurm["map_cluster"]["tasks-per-node"]), "--readFilesCommand", "gunzip", "-c", "--outSAMattributes", "All", "--outSAMreadID", "Number", "--outSAMtype", "BAM", "SortedByCoordinate", "--sjdbGTFfile", str(gene_gtf), "--genomeDir", str(star_index), "--limitBAMsortRAM", str(RAM)]
                if unique:
                    cmd.extend(["--outFilterMultimapNmax", "1", "--outFilterMismatchNoverLmax", "0.03"])
                else:
                    cmd.extend(["--outFilterMultimapNmax", "100", "--winAnchorMultimapNmax", "200"])
                if out_tmp_dir != None:
                    cmd.extend(["--outTmpDir", out_tmp_dir])
                
                cmd.extend(["--readFilesIn", os.path.join(fastq_dir, (self.cluster_name + "_R2.fastq.gz"))])
                if snic_tmp:
                    cmd.extend(["--outFileNamePrefix", (str(os.path.join("$SNIC_TMP/map_cluster/tmp", self.cluster_name)) + "_"), (" || exit 2; cp $SNIC_TMP/map_cluster/tmp/" + self.cluster_name + "_* " + outdir + "; echo 'Copying from tmp'")])
                else:
                    cmd.extend(["--outFileNamePrefix", (str(os.path.join(outdir, self.cluster_name)) + "_")])

                result = run_instruction(cmd = cmd, fun = "map_cluster", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = "map_cluster", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs["map_cluster"] = outdir
                return result
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)

    def gene_count(self, experiment_name, sample_id, bam, outdir, gene_gtf, s=1, unique=True, slurm=None, modules=None, snic_tmp = False, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("gene_count_scripts"):
                    os.makedirs("gene_count_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                
                if snic_tmp:
                    snic_bam_dir = "$SNIC_TMP/" # map_cluster/tmp/
                    cmd = ["cp", bam, os.path.join(snic_bam_dir, (self.cluster_name + "_Aligned.sortedByCoord.out.bam") + ";\n")]
                    bam = snic_bam_dir + self.cluster_name + "_Aligned.sortedByCoord.out.bam"
                else:
                    cmd = []

                if unique:
                    cmd.extend(["subread-parse.sh featureCounts", "-s", str(s), "-F", "GTF", "-g", "gene_name", "-a", gene_gtf, "-o", os.path.join(outdir, (experiment_name + "_" + self.cluster_name + "_uniqueMap.cntTable")), bam])
                    function_name = "gene_count"
                else:
                    msg = "ERROR: Gene counts available only for unique mapping configuration"
                    log.write(msg)
                    return 3

                result = run_instruction(cmd = cmd, fun = "gene_count", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = function_name, dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs[function_name] = outdir
                return result
                
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)


    def TE_count(self, experiment_name, sample_id, bam, outdir, gene_gtf, te_gtf, s=1, unique=False, slurm=None, modules=None, snic_tmp = False, dry_run = False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("TE_count_scripts"):
                    os.makedirs("TE_count_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
                
                if snic_tmp:
                    snic_bam_dir = "$SNIC_TMP/" # map_cluster/tmp/
                    cmd = ["cp", bam, os.path.join(snic_bam_dir, (self.cluster_name + "_Aligned.sortedByCoord.out.bam") + ";\n")]
                    bam = snic_bam_dir + self.cluster_name + "_Aligned.sortedByCoord.out.bam"
                else:
                    cmd = []

                if unique:
                    cmd.extend(["subread-parse.sh featureCounts", "-s", str(s), "-F", "GTF", "-g", "transcript_id", "-a", te_gtf, "-o", os.path.join(outdir, (experiment_name + "_" + self.cluster_name + "_uniqueMap.cntTable")), bam])
                else:
                    if s == 1:
                        stranded = "forward"
                    elif s == 2:
                        stranded = "reverse"
                    elif s == 0:
                        stranded = "no"
                        
                    cmd.extend(["TEcount", "-b", bam, "--GTF", gene_gtf, "--TE", te_gtf, "--format", "BAM", "--stranded", stranded, "--mode", "multi", "--sortByPos", "--project", os.path.join(outdir, (experiment_name + "_" + self.cluster_name + "_"))])
            
                if unique:
                    function_name = "TE_count_unique"
                else:
                    function_name = "TE_count"

                result = run_instruction(cmd = cmd, fun = "TE_count", name = ("sample_" + sample_id + "_cluster_" +  self.cluster_name), fun_module = function_name, dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules)
                exit_code = result[1]

                if exit_code == 0:
                    self.outdirs[function_name] = outdir
                return result
                
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                log.write(msg)


