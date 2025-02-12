#! /bin/python3
import subprocess
from subprocess import PIPE
import time
from .bcolors import Bcolors
import os

def run_instruction(cmd, fun, dry_run, fun_module, name, logfile, slurm = None, modules = None, modules_path = None):
    with open(logfile, "a") as log:
        log.write(' '.join(cmd))
        if slurm != None:
            cmd = ' '.join(cmd)
            job_file =  os.path.join((fun + "_scripts/"), (name + "_" + fun + ".sh"))
            try:
                if not dry_run:
                    job_id, exit_code, msg = run_job(fun_module, job_file, cmd, slurm, modules, modules_path)
                    # msg = sucess_submit(fun, name, job_id)
                    # log.write(msg)
                    
                    # exit_code = wait_for_job(job_id)
                    # msg = check_exit_codes(fun, name, job_id, exit_code)
                    log.write(msg)
                    log.write(str(job_id) + " returned " + str(exit_code))

                    return [msg, exit_code]
                else:
                    log.write(cmd + "\n")
                    run_job(fun_module, job_file, cmd, slurm, modules, dry_run)
                    return [cmd, 0]
            except:
                msg = generic_error(fun, name)
                log.write(msg)
                return [msg, 2]
        else:
            subprocess.call(cmd)

def run_job(function, job_file, code, slurm, modules, dry_run = False):
    with open(job_file, "w") as fout:
        fout.writelines("#!/bin/bash\n")
        try:
            for k,v in slurm[function].items():
                fout.writelines("#SBATCH --" + str(k) + "=" + str(v) + "\n")
        except:
            for k,v in slurm["__default__"].items():
                fout.writelines("#SBATCH --" + str(k) + "=" + str(v) + "\n")
        fout.writelines("module purge\n")
        try:
            if modules_path is not None:
                fout.writelines("module use " + modules_path + "\n")
            for i in modules[function]:
                fout.writelines("module load " + str(i) + "\n")
        except:
            pass
        fout.writelines('echo "' + code + '" || exit 2 \n')
        fout.writelines(code + " || exit 2 \n")

    if dry_run:
        return 

    sbatch_out = subprocess.run([("sbatch --wait " + str(job_file))], shell=True, stdout=PIPE, stderr=PIPE)
    time.sleep(3)
    job_id = (sbatch_out.stdout.split()[len(sbatch_out.stdout.split())-1]).decode("utf-8")
    exit_code = sbatch_out.returncode

    if exit_code == 0:
        msg = Bcolors.OKGREEN + job_id + " finished succesfully." + Bcolors.ENDC + "\n"
    elif exit_code == 1:
        msg = Bcolors.FAIL + " signal error." + Bcolors.ENDC + "\n"
    elif exit_code == 2:
        msg = Bcolors.FAIL + job_id + " finished with an error." + Bcolors.ENDC + "\n"
    return (job_id, exit_code, msg)

# def cancel(job_id):
#     scancel_out = subprocess.run([("scancel " + str(job_id))], shell=True, stdout=PIPE, stderr=PIPE)
#     time.sleep(3)
#     if scancel_out.returncode == 0:
#         print(Bcolors.FAIL + "Job" + str(job_id) + " got cancelled." + Bcolors.ENDC + "\n")
#         return
#     else:
#         print(Bcolors.FAIL + "Something went wrong. Please try again" + Bcolors.ENDC + "\n")
#         return

# def check_status(job_id):
#     job = subprocess.run(("sacct -j " + str(job_id) + " --format=state"), shell=True, stdout=PIPE, stderr=PIPE)
#     job = job.stdout.decode("utf-8").split()
#     status = job[(len(job)-1)]
#     return status.lower()

# def on_the_way(job_id):
#     status = check_status(job_id)
#     if status == "running" or status == "pending":
#         return True
#     else:
#         return False

# def wait_for_job(job_id):
#     # print("I'm here!")
#     time.sleep(3)
#     while on_the_way(job_id):
#         time.sleep(3)
#     if check_status(job_id) == "completed":
#         # print(Bcolors.OKGREEN + "Job " + job_id + " has been completed! :-)" + Bcolors.ENDC + "\n")
#         return 0
#     elif check_status(job_id) == "cancelled":
#         return 1
#     elif check_status(job_id) == "failed":
#         return 2
#     else:
#         return 3

# def check_exit_codes(fun, sample_id_cluster_name, job_id, exit_code):
#     if exit_code == 0:
#         return (Bcolors.OKGREEN + "Job " + str(job_id) + " finished " + fun + " succesfully. " + sample_id_cluster_name + Bcolors.ENDC + "\n")
#     elif exit_code == 1:
#         return (Bcolors.FAIL + "Job " + str(job_id) + " was cancelled during " + fun + ". " + sample_id_cluster_name + Bcolors.ENDC + "\n")
#     elif exit_code == 2:
#         return (Bcolors.FAIL + "Job " + str(job_id) + " failed during " + fun + ". " + sample_id_cluster_name + Bcolors.ENDC + "\n")
#     elif exit_code == 3:
#         return (Bcolors.FAIL + "Something strange happened to job " + str(job_id) + " during " + fun + ". " + sample_id_cluster_name + Bcolors.ENDC + "\n")

def generic_error(fun, info):
    return Bcolors.FAIL + "Something went wrong creating the " + fun + " job for " + info + Bcolors.ENDC + "\n"

# def sucess_submit(fun, info, job_id):
#     return Bcolors.OKBLUE + fun.capitalize() + " for " + info + " has been submitted and has job ID " + str(job_id) + Bcolors.ENDC + "\n"
