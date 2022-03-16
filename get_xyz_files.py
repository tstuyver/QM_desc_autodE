from argparse import ArgumentParser, Namespace
import tarfile
import os
import subprocess

parser = ArgumentParser()
parser.add_argument(
    "--tar_file", type=str, required=True, help="tar file to extract .xyz files from"
)

if __name__ == "__main__":
    args = parser.parse_args()
    if not os.path.isdir("xyz_files"):
        os.mkdir("xyz_files")
    path = os.path.join(os.getcwd(), "xyz_files")
    tar = tarfile.open(args.tar_file)
    os.chdir("/")
    tar.extractall()
    tar.close()
    os.chdir(path)
    dir_list = [path for path in os.listdir() if os.path.isdir(path)]
    for dir in dir_list:
        subprocess.run(["rm", "-r", f"{dir}"])
