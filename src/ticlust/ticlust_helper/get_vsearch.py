import platform
import os
import shutil
import requests
from pathlib import Path
import tarfile
import zipfile

global VSEARCH_REPO, BIN_DIR
VSEARCH_REPO = "torognes/vsearch"  # Replace with the actual repository
BIN_DIR = Path(__file__).parent.parent.joinpath("bin")
BIN_DIR.mkdir(parents=True, exist_ok=True)


def get_vsearch_latest() -> str:
    url = f"https://api.github.com/repos/{VSEARCH_REPO}/releases/latest"

    response = requests.get(url)
    if response.status_code == 200:
        latest_version = response.json()["tag_name"]
        return latest_version
    else:
        return "v2.29.3"  # Fallback to a known version


def get_vsearch_bin_url():
    system = platform.system().lower()
    arch = platform.machine().lower()
    latest_version = get_vsearch_latest()
    
    base_url = f"https://github.com/torognes/vsearch/releases/download/{latest_version}/"
    base_url = f"https://github.com/torognes/vsearch/releases/download/v2.29.3/"
    binaries = {
        "linux": {
            "aarch64": f"vsearch-{latest_version[1:]}-linux-aarch64.tar.gz",
            "x86_64": f"vsearch-{latest_version[1:]}-linux-x86_64.tar.gz",
            "riscv64": f"vsearch-{latest_version[1:]}-linux-riscv64.tar.gz",
            "mips64el": f"vsearch-{latest_version[1:]}-linux-mips64el.tar.gz",
            "ppc64le": f"vsearch-{latest_version[1:]}-linux-ppc64le.tar.gz"
        },
        "darwin": {
            "aarch64": f"vsearch-{latest_version[1:]}-macos-aarch64.tar.gz",
            "x86_64": f"vsearch-{latest_version[1:]}-macos-x86_64.tar.gz",
            "universal": f"vsearch-{latest_version[1:]}-macos-universal.tar.gz"
        },
        "windows": {
            "x86_64": f"vsearch-{latest_version[1:]}-win-x86_64.zip"
        }
    }
    
    if system in binaries and arch in binaries[system]:
        return base_url + binaries[system][arch]
    else:
        raise RuntimeError(f"Unsupported platform: {system} {arch}")


def download_vsearch_bin():
    binary_url = get_vsearch_bin_url()
    response = requests.get(binary_url)
    if response.status_code == 200:
        if 'text/html' in response.headers.get('Content-Type', ''):
            raise RuntimeError("Error: URL responded with HTML content instead of binary.")
        BIN_DIR.mkdir(parents=True, exist_ok=True)
        compressed_file = BIN_DIR / Path(binary_url).name
        with open(compressed_file, "wb") as f:
            f.write(response.content)
        return compressed_file
    else:
        return None


def extract_file_from(file_path, to_extract: str, extract_to="."):
    file_path = Path(file_path)
    extract_to = Path(extract_to)
    new_file = ''
    if file_path.suffix == ".zip":
        with zipfile.ZipFile(file_path, "r") as zip_ref:
            # extract only the file given in to_extract
            for file in zip_ref.namelist():
                if file.endswith(to_extract):
                    new_file = file_path.parent / file
                    zip_ref.extract(file, path=file_path.parent)
                    break
    elif file_path.suffix in [".tar", ".gz", ".tgz"]:
        with tarfile.open(file_path, "r:gz") as tar:
            for member in tar.getmembers():
                if member.name.endswith(to_extract):
                    # store the full path of file in the compressed file to extract it
                    new_file = file_path.parent / member.name
                    tar.extract(member, path=file_path.parent)
                    break
    else:
        raise RuntimeError(f"Unsupported archive format: {file_path}")
    # take the relative path of new_file to file_path.parent
    to_rem = new_file.relative_to(file_path.parent)
    # move the extracted file to the desired location
    if new_file:
        shutil.move(new_file, extract_to)
    else:
        raise RuntimeError(f"File {to_extract} not found in {file_path}")
    # remove to_rem
    shutil.rmtree(file_path.parent / to_rem.parts[0])
    return extract_to


def get_vsearch_bin_path() -> Path:
    vsearch_bin_zip = download_vsearch_bin()
    if not vsearch_bin_zip:
        raise RuntimeError("Failed to download the binary")    
    binary_path = vsearch_bin_zip.parent / "vsearch"
    extract_file_from(vsearch_bin_zip, "bin/vsearch", extract_to=binary_path)
    binary_path.chmod(0o775)
    vsearch_bin_zip.unlink()

    return binary_path


# Example usage:
if __name__ == "__main__":
    try:
        binary_path = get_vsearch_bin_path()
        print(f"Binary downloaded to: {binary_path}")
    except RuntimeError as e:
        print(e)