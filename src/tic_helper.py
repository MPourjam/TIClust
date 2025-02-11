# pylint: disable=missing-docstring,missing-class-docstring,missing-function-docstring
# pylint: disable=too-few-public-methods,no-method-argument,too-many-lines
# pylint: disable=logging-fstring-interpolation,too-many-arguments,too-many-locals
import re
import logging
import subprocess
import gzip
import shutil
import zipfile
import inspect
import threading
import tempfile
from pathlib import Path, PurePath
from datetime import datetime as dt
import mimetypes as mtypes
from os import (
    path as ospath,
    remove as os_remove,
    replace,
    access,
    X_OK
)
try:
    # Posix based file locking (Linux, Ubuntu, MacOS, etc.)
    #   Only allows locking on writable files, might cause
    #   strange results for reading.
    import fcntl

    def lock_file(f):
        if f.writable():
            fcntl.lockf(f, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def unlock_file(f):
        if f.writable():
            fcntl.lockf(f, fcntl.LOCK_UN)
except ModuleNotFoundError:
    # Windows file locking
    import msvcrt

    def file_size(f):
        return ospath.getsize(ospath.realpath(f.name))

    def lock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_RLCK, file_size(f))

    def unlock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, file_size(f))


class ArgsetException(Exception):
    pass


# overwriting system_sub() to run it with subprocess.run
def loud_subprocess(cmd_args_list: list, shell_bool: bool = False, cap_output: bool = False):
    dev_null_rgx = re.compile(r'(\s+)?([12]?>)(\s+)?(/dev/null|&1|&2)')
    if not isinstance(cmd_args_list, list):
        raise ValueError("cmd_args_list should be a list")
    cmd_string = "\t\t".join(cmd_args_list)
    # remove any null redirector from the given command
    # remove  > /dev/null 2>&1
    cmd_string = re.sub(dev_null_rgx, "", cmd_string)
    cmd_list = cmd_string.split("\t\t")
    # If capture_output_bool is true then do not redirect to /dev/null
    run_output = subprocess.run(
        cmd_list,
        capture_output=cap_output,
        encoding='utf-8',
        shell=shell_bool,
        check=True
    )

    return run_output, cmd_list


class LogPrint(logging.Logger):

    def __init__(self, name):
        super().__init__(name)

    def debug(self, msg, *args, **kwargs):
        print(f"\033[95m{msg}\033[0m")

    def info(self, msg, *args, **kwargs):
        print(f"\033[94m{msg}\033[0m")

    def warning(self, msg, *args, **kwargs):
        print(f"\033[93m{msg}\033[0m")

    def error(self, msg, *args, **kwargs):
        print(f"\033[91m{msg}\033[0m")

    def critical(self, msg, *args, **kwargs):
        print(f"\033[1m\033[91m{msg}\033[0m")


def system_sub(
    cmd_args_list: list,
    force_log: bool = False,
    shell: bool = False,
    capture_output: bool = True,
    quiet: bool = False,
    logger_obj: logging.Logger = LogPrint("print_logger")):
    # logging
    if force_log:
        msg = f"COMMAND: {' '.join(cmd_args_list)}\n\n"
        logger_obj.debug(msg)

    run_output, cmd_list = loud_subprocess(
        cmd_args_list,
        shell_bool=shell,
        cap_output=capture_output
    )
    if run_output.returncode == 137:  # Process killed due to memory limit
        err_msg = f"Command {' '.join(cmd_list)} exceeded memory limit."
        err_msg += "\n\tIf you are using usearch 32-bit version, consider"\
                    " upgrading to 64-bit version."
        err_msg += "\n\tIf you are using usearch 64-bit then run the programm"\
                    " with lower number of threads."
        raise MemoryError(err_msg)
    if run_output.returncode != 0 and not quiet:
        raise Exception(
            f"Error in running command {' '.join(cmd_list)}:\n{run_output.stderr}"
        )
    if run_output.returncode != 0 and quiet:
        logger_obj.error(
            f"Error in running command {' '.join(cmd_list)}:\n{run_output.stderr}"
        )
    elif run_output.returncode == 0 and run_output.stderr and not quiet:
        logger_obj.debug(
            f"Warning in running command {' '.join(cmd_list)}:\n{run_output.stderr}"
        )
    elif run_output.returncode == 0 and not quiet:
        logger_obj.debug(
            f"Logs in running command {' '.join(cmd_list)}:\n{run_output.stdout}"
        )
    return run_output


def gzip_to_fastq(*files) -> list:
    '''
    It takes files gunzip or zip them and returns the name with proper fastq fuffix.
    '''
    file_path = [Path(PurePath(f)) for f in files if f]
    file_path = [f for f in file_path if f.is_file()]
    if len(files) != len(file_path):
        raise TypeError("Some given arguments are not files.")
    fastq_paths = []
    for f in file_path:
        app, typ = mtypes.guess_type(f)
        file_dir, file_name = f.parent, str(f.name)
        file_name = file_name.replace(".", "_").replace("_gz", "").replace("_bz2", "")
        asciifile_name = file_name.replace(" ", "").replace("_fastq", "") + ".fastq"
        asciifile_name_normalized = ''.join(
            e
            for e in asciifile_name
            if e.isalnum() or e in ['_', '-', '.']
        )
        asciifile = str(file_dir.joinpath(asciifile_name_normalized))
        # derefrencing f if it's a link
        f_real = f.resolve()
        if 'zip' in str(typ):
            # For gzipped files
            with gzip.open(f_real, 'rb') as gz_file:
                with open(asciifile, 'wb') as ascii_file:
                    shutil.copyfileobj(gz_file, ascii_file)
                # Remove the original gzip file
            os_remove(f)
            fastq_paths.append(asciifile)
        elif 'zip' in str(app):  # For zipped files
            with zipfile.ZipFile(f_real, 'r') as zip_ref:
                zip_ref.extractall(asciifile)
            # Remove the original zip file
            os_remove(f)
            fastq_paths.append(asciifile)
        else:
            try:
                with open(f_real, "r+", encoding="utf-8") as fi:
                    line = fi.readline()
                if len(line) == len(line.encode()):  # If it's ascii
                    try:
                        shutil.copy(f_real, asciifile)
                    except shutil.SameFileError:
                        pass
                    fastq_paths.append(asciifile)
            except Exception as e:
                print(f"{f}\t{e}")
                raise e
    if len(file_path) != len(fastq_paths):
        raise FileNotFoundError("Some files could not get converted or were not in utf-8 format!")
    # returns absolute paths
    return fastq_paths


def find_files_and_dirs_owned_by_root(directory):
    root_files_and_dirs = []
    directory = Path(directory).absolute()

    for path in directory.rglob('*'):
        try:
            # Get file or directory owner information
            if path.is_symlink():
                raise TypeError("Symlinks are skipped.")
            file_stat = path.stat()
            # Check if the owner is root (uid 0)
            if file_stat.st_uid == 0:
                root_files_and_dirs.append(str(path))
        except TypeError as texc:
            print(f"Error while checking for permissions on {path}: {texc}")
        except Exception:
            pass

    return root_files_and_dirs


def generate_timestamp(thread_safe=False):
    current_time = dt.now()
    timestamp = current_time.strftime("%Y%m%d_%H%M%S")
    if thread_safe:
        # Get the current thread obj hex
        thread_obj = threading.current_thread()
        thread_obj = str(hash(id(thread_obj)))
        timestamp += f"_{thread_obj}"

    return timestamp


def get_file_handler(logger, file_path) -> logging.FileHandler:
    abs_file_path = Path(file_path).resolve()
    fh = logging.FileHandler(abs_file_path, mode='w')
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            if Path(handler.baseFilename).resolve() == file_path:
                return handler
    return fh


def get_stream_handler(logger) -> logging.StreamHandler:
    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler):
            return handler
    return logging.StreamHandler()


def gimmelogger(
    logger_name: str = "",
    log_file: str = "",
    only_file: bool = True,
    propagate: bool = True
    ):
    # finding caller file name and setting logger file path
    caller_frame = inspect.currentframe().f_back
    logger_name = Path(
        PurePath(caller_frame.f_code.co_filename)
        ).stem if not logger_name else str(logger_name)
    if not log_file:
        log_file_path = Path(
            PurePath(
                inspect.getframeinfo(caller_frame).filename
                )
            ).parent.joinpath(f"{logger_name}_log.txt")
    else:
        log_file_path = Path(PurePath(log_file)).absolute()

    # Setting logger formatter with line number of source of log
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s:%(filename)s:%(lineno)d - %(levelname)s - %(message)s'
    )

    logger = logging.getLogger(logger_name)
    # set the logging level
    logger.setLevel(logging.DEBUG)
    if not only_file:
        ch = get_stream_handler(logger)
        ch.setLevel(logging.INFO)
        # create a formatter
        # set the formatter to the console handler
        ch.setFormatter(formatter)
        # add the console handler to the logger
        logger.addHandler(ch)
    # Logger
    if log_file:
        if not log_file_path.parent.is_dir():
            log_file_path.parent.mkdir(parents=True, exist_ok=True)
        fh = get_file_handler(logger, log_file_path)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    # Propagate option
    logger.propagate = propagate
    return logger


def slice_list(li: list, n: int) -> list:
    return [li[i:i + n] for i in range(0, len(li), n)]


class FileUtil:

    @staticmethod
    def is_gzip(file_path):
        file_path = Path(PurePath(file_path)).absolute()
        try:
            with gzip.open(file_path, 'rb') as f:
                f.read(1)
            return True
        except OSError:
            return False

    @staticmethod
    def is_zip(file_path):
        file_path = Path(PurePath(file_path)).absolute()
        try:
            with zipfile.ZipFile(file_path, 'r') as f:
                f.testzip()
            return True
        except zipfile.BadZipFile:
            return False

    @staticmethod
    def gunzip(file_path: str) -> str:
        """
        It gunzips the file and returns the path of the binary file.
        """
        file_path = Path(PurePath(file_path)).absolute()
        with gzip.open(file_path, 'rb') as f:
            file_content = f.read()
        bin_file = Path(PurePath(file_path).parent).joinpath(Path(PurePath(file_path).stem))
        with open(bin_file, 'wb') as f:
            f.write(file_content)
        if FileUtil.is_binary(bin_file):
            return str(bin_file)
        return ""

    @staticmethod
    def is_binary(file_path: str) -> bool:
        file_path = Path(PurePath(file_path)).absolute()
        if not file_path.is_file():
            return False
        # Check if the file is a binary file
        with open(file_path, "rb") as f:
            header = f.read(4)
            if header != b'\x7fELF':
                return False
        return True

    @staticmethod
    def change_mode(file_path: str, mode: int = 0o777):
        file_path = Path(PurePath(file_path)).absolute()
        file_path.chmod(mode)
        return mode

    @staticmethod
    def is_executable(file_path: str) -> bool:
        file_path = Path(PurePath(file_path)).absolute()
        if not file_path.is_file() or not access(file_path, X_OK):
            return False
        return True


class Usearch(FileUtil):

    desired_version = "usearch v11"

    def __init__(self, file: str, version: str = ""):
        self.file = Path(PurePath(file)).absolute()
        self.binfile = ""
        if not self.file.is_file():
            raise FileNotFoundError(f"File {file} does not exist.")
        self.desired_version = version if version else self.desired_version

    def which_usearch_version(self) -> str:
        """
        It returns the version of the usearch binary file.
        """
        version_cmd = [str(self.binfile), "--version"]
        version_out = subprocess.run(
            version_cmd,
            capture_output=True,
            text=True,
            check=True
        )
        return version_out.stdout

    def is_version(self, version: str = "") -> bool:
        """
        It checks if the given version is the same as the desired version.
        """
        version = version if version else self.desired_version
        return version in self.which_usearch_version()

    def check_or_get_bin(self) -> tuple:
        """
        It checks the given file path and if it's a gzip file then it looks for a linux binary file
         made for i86 architecture and returns the path of the binary file.
        In case the file is not a gzip then it checks if the file is binary and made for i86
         architecture and returns the path of the binary file.
        """
        ret_tup = (1, "")
        if FileUtil.is_binary(self.file):
            # if the file is a binary file
            self.binfile = self.file
        elif FileUtil.is_gzip(self.file):
            # If the file is a gzip file
            self.binfile = FileUtil.gunzip(self.file)
        else:
            # If the file is not a binary file
            # then we check if the binary file is inside
            bin_file = self.file.parent.joinpath(self.file.stem)
            if FileUtil.is_binary(bin_file):
                self.binfile = bin_file
            else:
                # if the binary file is not found
                self.binfile = ""

        if self.binfile:
            # if the binary file is found then we make it executable
            FileUtil.change_mode(self.binfile, 0o777)
            # if the version is not the desired version then we return 162
            ret_tup = (162, "")
            if self.is_version():
                ret_tup = (0, self.binfile)
        else:
            # if the binary file is not found then we return 161
            ret_tup = (161, "")

        return ret_tup

    def __bool__(self):
        return self.binfile != ""

    def __str__(self):
        return str(self.file)

    def __repr__(self):
        return f"Usearch {self.desired_version}" if self.binfile else "Usearch"

def onelinefasta(fastafilepath):
    proper_filepath = Path(fastafilepath).resolve()
    dirpath = proper_filepath.parent
    filename = proper_filepath.name
    with tempfile.NamedTemporaryFile(
        delete=False,
        suffix=str(proper_filepath.suffix),
        dir=str(dirpath)
    ) as temp_file:
        newfile_temp_name = temp_file.name
    with proper_filepath.open(
            'r',
            encoding='utf-8'
        ) as f, open(
            newfile_temp_name,
            "w+", encoding='utf-8'
        ) as onelinefa:
        line = f.readline()
        sequence = ""
        while line:
            if line[0] == '>':
                if sequence:
                    onelinefa.write(sequence + "\n")
                    sequence = ""
                onelinefa.write(line)
            elif line in ['\n', '\r\n']:
                pass
            else:
                sequence += line.strip()
            line = f.readline()
        if sequence:
            onelinefa.write(sequence + "\n")
    try:
        replace(newfile_temp_name, filename)
    finally:
        if Path(newfile_temp_name).exists():
            Path(newfile_temp_name).unlink()
