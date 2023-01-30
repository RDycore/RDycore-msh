
import os
import zipfile
from github_release import gh_asset_download

HERE = os.path.abspath(os.path.dirname(__file__))

if (__name__ == "__main__"):
    """
    Download and unzip assets from rdycore-msh release.

    """
    cwd_pointer = os.getcwd()

    try:
        datadir = os.path.join(HERE, "data")
        os.chdir(datadir)

        print("Extracting to:", datadir)
        
        repo = "RDycore/RDycore-msh"
        vtag = "v0.0.1"
        gh_asset_download(repo, vtag)

        for item in os.listdir(datadir):
            if item.endswith(".zip"):
                filename = os.path.abspath(item)
                zip_func = zipfile.ZipFile(filename)
                zip_func.extractall(datadir)
                zip_func.close()

        print("")
        print("Done downloading + unpacking datasets.")
        print("")

    finally:
        os.chdir(cwd_pointer)



