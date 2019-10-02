import os
import tarfile
import urllib.request
import zipfile
import geopandas as gpd

def download_to_path(path: str, url: str):
    r"""
    Download from a HTTP or FTP url to a filepath.

    >>> d = download_to_path(
    ...    path="data/Arthern_accumulation_bedmap2_grid1.tif",
    ...    url="https://secure.antarctica.ac.uk/data/bedmap2/resources/Arthern_accumulation/Arthern_accumulation_tif.zip",
    ... )
    """

    folder, filename = os.path.split(p=path)
    downloaded_filename = os.path.basename(urllib.parse.urlparse(url=url).path)

    # Download file using URL first
    if not os.path.exists(os.path.join(folder, downloaded_filename)):
        r = urllib.request.urlretrieve(
            url=url, filename=os.path.join(folder, downloaded_filename)
        )

    # If downloaded file is not the final file (e.g. file is in an archive),
    # then extract the file from the archive!
    if filename != downloaded_filename:
        # Extract tar.gz archive file
        if downloaded_filename.endswith(("tgz", "tar.gz")):
            try:
                archive = tarfile.open(name=f"{folder}/{downloaded_filename}")
                archive.extract(member=filename, path=folder)
            except:
                raise
        # Extract from .zip archive file
        elif downloaded_filename.endswith((".zip")):
            try:
                archive = zipfile.ZipFile(file=f"{folder}/{downloaded_filename}")
                archive.extract(member=filename, path=folder)
            except:
                raise
        else:
            raise ValueError(
                f"Unsupported archive format for downloaded file: {downloaded_filename}"
            )

    return os.path.exists(path=path)

