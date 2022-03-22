from argparse import ArgumentParser, Namespace
import tarfile
import os
import subprocess

parser = ArgumentParser()
parser.add_argument(
    "--tar_file", type=str, required=True, help="tar file to extract .xyz files from"
)
parser.add_argument("--folder_name", type=str, default="XTB_opt", help="name of the folder to which the tar-archive is extracted")

if __name__ == "__main__":
    args = parser.parse_args()
    #if not os.path.isdir("xyz_files"):
    #    os.mkdir("xyz_files")
    path = os.path.join(os.getcwd(), "xyz_files")
    path_tar = os.path.join(os.getcwd(), args.tar_file)
    path_folder = os.path.join(os.getcwd(), args.folder_name)
    tar = tarfile.open(path_tar)
    os.chdir("/")
    tar.extractall()
    tar.close()
    #os.chdir(path_folder)
    dir_list = [item for item in os.listdir(path_folder) if os.path.isdir(item)]
    for dir in dir_list:
        subprocess.run(["rm", "-r", f"{dir}"])
    subprocess.run(["mv", "-r", path_tar, path])
