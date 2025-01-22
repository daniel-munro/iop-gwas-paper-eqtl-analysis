import os
import pysftp

os.chdir("/gpfs/home/dmunro/data/eyes/fastq")

files = open("/gpfs/home/dmunro/data/eyes/todo.txt", "r").read().splitlines()

with pysftp.Connection("transfer.uthsc.edu", username="redacted", password="redacted") as sftp:
    with sftp.cd("/HS_eyes"):
        # sftp.get_r("C202SC19070808_1", ".")
        # sftp.get_r("C202SC19070808_2", ".")
        # sftp.get_r("HS_eyes/RSEM_outputs", ".")
        for fname in files:
            print(fname)
            # dr = fname.split("/")[0]
            sftp.get(fname, fname)
