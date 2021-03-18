from sys import argv
from ftplib import FTP
import gzip
import shutil
import os

##Collect command-line options in a dictionary
def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0][1:]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

args_oi = getopts(argv)

##These are file types that can downloaded and resp extensions
file_suffix_types_dict = dict()
file_suffix_types_dict["gb_fn"] = "_genomic.gbff.gz"
file_suffix_types_dict["ass_fn"] = "_assembly_report.txt"
file_suffix_types_dict["fna_fn"] = "_genomic.fna.gz"
file_suffix_types_dict["ft_fn"] = "_feature_table.txt.gz"
file_suffix_types_dict["cds_fn"] = "_cds_from_genomic.fna.gz" 
file_suffix_types_dict["prot_fn"] = "_protein.faa.gz"
file_suffix_types_dict["rna_fn"] = "_rna_from_genomic.fna.gz"
file_suffix_types_dict["ftcnt_fn"] = "_feature_count.txt.gz"

##A list of available files
available_files=[fn_oi for fn_oi in os.listdir(args_oi["od"]) if ".gz" not in fn_oi]

##Login for NCBI FTP
ftp = FTP('ftp.ncbi.nih.gov')
ftp.login('anonymous', 'saurabh.mk@gmail.com') ##change abc@abc.com to email id


downloaded_fna_list=[]
not_found_list=[]

with open(args_oi["wd"]+args_oi["list_fn"], "r") as required_genomes_details_fh:
    ##iterate over required assembly accessions ids
    for line in required_genomes_details_fh.readlines():
        ##only if accession contains string GCF/GCA
        if((line.strip()!="NA") & (line.strip()[0:3] in ["GCF", "GCA"])):
            
            ftp_folder_strings = line.strip().split("_")[0:2]
            ##Create a string for required genome's folder at FTP           
            ftp_folder = ftp_folder_strings[0]+"/"+ftp_folder_strings[1][0:3]+"/"+ftp_folder_strings[1][3:6]+"/"+ftp_folder_strings[1][6:9]+"/"
            ##change directory on FTP to that folder
            ftp.cwd('/genomes/all/'+ftp_folder)
            ##get a list of all versions available
            av_versions = ftp.nlst()
            ##get the latest version which should be listed last
            ftp.cwd(av_versions[0])
            GCF_id=av_versions[0]
            ##if this file is already available
            if(any(GCF_id in file_oi for file_oi in available_files)):
                with open(args_oi["sd"]+args_oi["fn_type"]+"_GCFids_found.txt", "a") as f_fh: #log download in the clade specific list
                    f_fh.write("%s\n" % GCF_id)
                print("you have the file already!")
                continue ##go to next file via for loop
            
            ##create string fo rexact file required with extension
            fn_oi = GCF_id+file_suffix_types_dict[args_oi["fn_type"]]

            print("looking for files in: " + '/genomes/all/' + ftp_folder + av_versions[0])
            try:
                ftp.cwd('/genomes/all/'+ftp_folder+GCF_id+"/")
                #download the file
                with open(args_oi["od"]+fn_oi, 'wb') as fh_oi:
                    ftp.retrbinary('RETR %s' % fn_oi, fh_oi.write)
                #write accession in list of files found
                with open(args_oi["sd"]+args_oi["fn_type"]+"_GCFids_found.txt", "a") as f_fh: #log download in the clade specific list
                    f_fh.write("%s\n" % GCF_id)
                
                #as long as file is not assembly summary, we need to decompress it
                if(args_oi["fn_type"] not in ["ass_fn"]):
                    with gzip.open(args_oi["od"]+fn_oi, "rb") as zh_in:
                        with open(args_oi["od"]+fn_oi[0:-3], "wb") as zh_out:
                            shutil.copyfileobj(zh_in, zh_out)
                    os.remove(args_oi["od"]+fn_oi)
            ##if file is not found
            except Exception as e:
                with open(args_oi["sd"]+args_oi["fn_type"]+"_GCFids_found.txt", "a") as f_fh:
                    f_fh.write("NA\n")
                print(e)
        else:
            print("not found")
            with open(args_oi["sd"]+args_oi["fn_type"]+"_GCFids_found.txt", "a") as f_fh:
                f_fh.write("NA\n")
            
ftp.quit()
