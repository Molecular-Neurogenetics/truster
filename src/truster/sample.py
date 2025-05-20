#! /bin/python3
import subprocess
import os
import json
import re
import sys
from .jobHandler import *
import pysam
import concurrent.futures
import pandas as pd
from .cluster import Cluster
from .bcolors import Bcolors

class Sample:
    def __init__(self, slurm, modules, logfile, sample_id="", sample_name="", raw_path = "", rdata_path = ""):
        self.slurm = slurm
        self.modules = modules
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.raw_path = raw_path
        self.rdata_path = rdata_path
        self.logfile = logfile

    def md5sum(self, outdir):
        if not os.path.exists("md5sum_scripts/"):
            os.makedirs("md5sum_scripts", exist_ok=True)
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)

        samples_files = []
        sample_outfile = os.path.join(outdir, (self.sample_name + "_md5.txt"))
        with open(self.logfile, "a") as log:
            try:
                for root, dirs, files in os.walk(self.raw_path):
                  for file in files:
                    sample_name_format = self.sample_name.replace("-", "_")
                    log.write("Checking for files: " + self.sample_name + " or, " + sample_name_format + "\n")
                    if self.sample_name in file or sample_name_format in file:
                       samples_files.append(file)
                       cmd = "md5sum " + os.path.join(root, file) + " >> " + sample_outfile + " || exit 2"
                       subprocess.run(cmd, shell=True)
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def quantify(self, cr_index, indir, outdir, nuclei = False):
        dry_run = False
        if not os.path.exists("quantify_scripts/"):
            os.makedirs("quantify_scripts", exist_ok=True)
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        sample_outdir = os.path.join(outdir, self.sample_id)
        with open(self.logfile, "a") as log:
            try:
                cmd = ["cellranger count"]
                if nuclei:
                    cmd.append("--include-introns true")
                cmd.extend(["--id", self.sample_id, "--transcriptome", cr_index, "--create-bam true --sample", self.sample_id, "--fastqs", indir, " || exit 2\nmv", self.sample_id, sample_outdir])
                log.write("About to run instruction")
                result = run_instruction(cmd = cmd, fun = "quantify", fun_module = "quantify", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]
                if exit_code == 0:
                    self.quantify_outdir = sample_outdir
                return exit_code
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def set_quantification_outdir(self, cellranger_outdir):
        self.quantify_outdir = cellranger_outdir
        with open(self.logfile, "a") as log:
            msg = "Quantification directory for sample" + self.sample_id + " is set to: " + cellranger_outdir + ".\n"
            log.write(msg)

    def velocity(self, te_gtf, gene_gtf, indir, dry_run=False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("velocity_scripts/"):
                    os.makedirs("velocity_scripts", exist_ok=True)
                
                cmd = ["velocyto run10x", "-m", te_gtf, indir, gene_gtf]
                result = run_instruction(cmd = cmd, fun = "velocity", fun_module = "velocity", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]
                
                if exit_code == 0:
                    subprocess.call("mv", self.sample_id, indir)
                
                return exit_code
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def plot_velocity(self, loom, indir, outdir):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("velocity_scripts/"):
                    os.makedirs("velocity_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)

                cwd = os.path.dirname(os.path.realpath(__file__))
                os.path.join(cwd, "py_scripts/plot_velocity.py")
                # print('plot_velocity -l <loom> -n <sample_name> -u <umap> -c <clusters> -o <outdir>')
                cmd = ["python", os.path.join(cwd, "py_scripts/plot_velocity"), "-l", loom, "-n", self.sample_id, "-u", os.path.join(indir, (self.sample_id + "_cell_embeddings.csv")), "-c", os.path.join(indir, (self.sample_id + "_clusters.csv")), "-o", outdir]
                result = run_instruction(cmd = cmd, fun = "plot_velocity", fun_module = "plot_velocity", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]
                return exit_code
            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)

    def empty_clusters(self):
        self.clusters = []

    def get_clusters(self, outdir, res = 0.5, perc_mitochondrial = None, min_genes = None, max_genes = None, normalization_method = "LogNormalize", max_size=500, dry_run = False):
        
        with open(self.logfile, "a") as log:
            if hasattr(self, "quantify_outdir"):
                try:
                    msg = f"Running get_clusters() for sample {self.sample_id}\n"
                    log.write(msg)
                    print(msg.strip())

                    if not os.path.exists("get_clusters_scripts"):
                        os.makedirs("get_clusters_scripts", exist_ok=True)
                    if not os.path.exists(outdir):
                        os.makedirs(outdir, exist_ok=True)

                    res = str(res)
                    max_size = str(max_size)
                    
                    cwd = os.path.dirname(os.path.realpath(__file__))
                    cmd = [os.path.join(cwd, "r_scripts/get_clusters.R"), "-i", os.path.join(self.quantify_outdir, "outs/filtered_feature_bc_matrix"), "-o", outdir, "-s", self.sample_id, "-r", res, "-n", normalization_method, "-S", max_size]

                    if perc_mitochondrial is not None:
                        cmd.extend(["-p", str(perc_mitochondrial)])
                    if min_genes is not None:
                        cmd.extend(["-m", str(min_genes)])
                    if max_genes is not None:
                        cmd.extend(["-M", str(max_genes)])

                    result = run_instruction(cmd = cmd, fun = "get_clusters", fun_module = "get_clusters", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                    exit_code = result[1]

                    if exit_code == 0:
                        self.clusters = [Cluster(j.split(".tsv")[0], os.path.join(outdir, j), self.logfile) for j in os.listdir(outdir) if j.endswith(".tsv")]
                        self.rdata_path = os.path.join(outdir, (self.sample_id + ".rds"))
                    return exit_code

                except KeyboardInterrupt:
                    msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                    log.write(msg)
            else:
                msg = f"get_clusters() error: Quantification directory not found for sample {self.sample_id}\n"
                log.write(msg)
                print(f"{Bcolors.FAIL}{msg.strip()}{Bcolors.ENDC}")
                return 4

    def register_clusters_from_path(self, path):
        self.clusters = [Cluster(j.split(".tsv")[0], os.path.join(path, j), self.logfile) for j in os.listdir(path) if j.endswith(".tsv")]
        self.rdata_path = os.path.join(path, (self.sample_id + ".rds"))
        msg = "Clusters from " + self.sample_id + " have been registered from " + path + ". Clusters: " + ', '.join([i.cluster_name for i in self.clusters]) + "\n"
        return msg

    def normalize_TE_counts(self, indir, outdir, dry_run=False):
        with open(self.logfile, "a") as log:
            try:
                if not os.path.exists("normalize_TE_counts_scripts"):
                    os.makedirs("normalize_TE_counts_scripts", exist_ok=True)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, exist_ok=True)
    
                cwd = os.path.dirname(os.path.realpath(__file__))
                cmd = ["Rscript", os.path.join(cwd, "r_scripts/normalize_TEexpression.R"), "-m", "individual", "-o", outdir, "-i", indir, "-r", self.rdata_path, "-n", self.sample_id]
                
                result = run_instruction(cmd = cmd, fun = "normalize_TE_expression", fun_module = "normalize_TE_expression", dry_run = dry_run, name = self.sample_id, logfile = self.logfile, slurm = self.slurm, modules = self.modules)
                exit_code = result[1]
                        
                if exit_code == 0:
                    self.norm_TE_counts = os.path.join(outdir, "TE_normalizedValues_aggregatedByClusters_melted.csv")
                return exit_code

            except KeyboardInterrupt:
                msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC + "\n"
                log.write(msg)


    def tsv_to_bam_clusters(self, outdir, slurm=None, modules=None, modules_path = None, remote_bam_flag = False, remote_bam_paths = None, remote_bam_tmpdir = None, dry_run = False):
         with open(self.logfile, "a") as log:
            outdir_sample = os.path.join(outdir, "tsv_to_bam/", (self.sample_id + "/"))
            if not os.path.exists(outdir_sample):
                os.makedirs(outdir_sample, exist_ok=True)
            
            if not remote_bam_flag:
                if os.path.isdir(self.quantify_outdir):
                    bam = os.path.join(self.quantify_outdir, "outs/possorted_genome_bam.bam")
                else:
                    msg = "Error: File not found. Please make sure that " + sample.quantify_outdir + " exists.\n"
                    log.write(msg)
                    return 4
                try: 
                    tsvs = ' '.join([c.tsv for c in self.clusters])
                    cmd = ["Rustody multi_subset_bam", "-t", "CB", "--bam", bam, "--values", tsvs, "--ofile", outdir_sample]
                    log.write(" ".join(cmd) + "\n\n")
                    result = run_instruction(cmd = cmd, fun = "tsv_to_bam_clusters", name = ("sample_" + self.sample_id), fun_module = "tsv_to_bam_clusters", dry_run = dry_run, logfile = self.logfile, slurm = slurm, modules = modules, modules_path = modules_path)
                    return result
                except KeyboardInterrupt:
                    msg = Bcolors.HEADER + "User interrupted" + Bcolors.ENDC
                    log.write(msg)
            # We have a remote bam file
            # But I dont want to load to the node the bam files one per cluster in a sample, so, I'll put them all in one job
            # For my use case this is ok (dealing with 5 clusters), but beware!
            else: 
                remote_bam_tmpdir = os.path.join(remote_bam_tmpdir, (self.sample_id + ".bam"))
                bam_paths = pd.read_csv(remote_bam_paths, sep = "\t")
                bam_paths = bam_paths.set_index("sample_id")
                sample_remote_bam_path = bam_paths.loc[self.sample_id, "path_bam"]
                bam = remote_bam_tmpdir # We will move it to the tmp dir, so we can write the command directly pointing there
                bai = (remote_bam_tmpdir + ".bai")

                # Let's copy it to its tmp dir using iget (bam path was already fixed in experiment.py tsv_to_bam_clusters()) 
                cmd = ["iget", sample_remote_bam_path, bam, " || exit 2; ", "iget", (sample_remote_bam_path + ".bai"), bai, " || exit 2; "] 
                
                for cluster in self.clusters:
                    cluster_cmd = ["subset-bam", "--bam", bam, "--cell-barcodes", cluster.tsv, "--out-bam", os.path.join(outdir_sample, (cluster.cluster_name + ".bam || echo FAIL: " + cluster.cluster_name + "; "))]    
                    cmd.extend(cluster_cmd)
                log.write(" ".join(cmd) + "\n\n")
                result = run_instruction(cmd = cmd, fun = "tsv_to_bam_clusters", dry_run = dry_run, fun_module = "tsv_to_bam_clusters", name = ("sample_" + self.sample_id), logfile = self.logfile, slurm = slurm, modules = modules, modules_path = modules_path)
                return result