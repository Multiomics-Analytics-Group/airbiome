import logging
import os
import sys
import yaml
from datetime import datetime
from pathlib import Path


def load_config(config_path):  # <-- load_config is defined here
    """
    Loads configuration from a YAML file.
    ...
    """
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# --- LOGGING UTILITY FUNCTIONS ---
def get_basename(fname: None | str = None) -> str:
    """
    - For a given filename, returns basename WITHOUT file extension
    - If no fname given (i.e., None) then return basename that the function is called in

    Parameters
    ----------
    fname: str, optional
        The filename to get basename from. Default is None.

    Returns
    -------
    str
        basename of given filepath or the current file the function is executed

    Examples
    ---------
    1)
    >>> get_basename()
    utils

    2)
    >>> get_basename('this/is-a-filepath.csv')
    is-a-filepath
    """
    if fname is not None:
        # PRECONDITION
        # Assuming check_path is available or replaced by Path(fname).exists()
        if not Path(fname).exists():  # Replaced check_path with pathlib's exists()
            raise FileNotFoundError(f"The specified path does not exist: {fname}")
        # MAIN FUNCTIONS
        return os.path.splitext(os.path.basename(fname))[0]
    else:
        return os.path.splitext(os.path.basename(sys.argv[0]))[0]


def get_time(incl_time: bool = True, incl_timezone: bool = True) -> str:
    """
    Gets current date, time (optional) and timezone (optional) for file naming

    Parameters
    ----------
    - incl_time (bool): whether to include timestamp in the string
    - incl_timezone (bool): whether to include the timezone in the string

    Returns
    -------
    str
        fname that includes date, timestamp and/or timezone
        connected by '_' in one string e.g.YYYYMMdd_hhmm_timezone

    Examples
    --------
    1)
    >>> get_time()
    '20231019_101758_CEST'

    2)
    >>> get_time(incl_time=False)
    '20231019_CEST'

    """

    # PRECONDITIONALS
    assert isinstance(incl_time, bool), "incl_time must be True or False"
    assert isinstance(incl_timezone, bool), "incl_timezone must be True or False"

    # MAIN FUNCTION
    # getting current time and timezone
    the_time = datetime.now()
    timezone = datetime.now().astimezone().tzname()
    # convert date parts to string

    # putting date parts into one string
    if incl_time and incl_timezone:
        fname = the_time.isoformat(sep="_", timespec="seconds") + "_" + timezone
    elif incl_time:
        fname = the_time.isoformat(sep="_", timespec="seconds")
    elif incl_timezone:
        fname = "_".join([the_time.isoformat(sep="_", timespec="hours")[:-3], timezone])
    else:
        y = str(the_time.year)
        m = str(the_time.month)
        d = str(the_time.day)
        fname = y + m + d

    # optional
    fname = fname.replace(":", "-")  # remove ':' from hours, minutes, seconds

    return fname


def generate_log_filename(folder: str = "logs", suffix: str = "") -> str:
    """
    Creates log file name and path

    Parameters
    ----------
    folder (str): name of the folder to put the log file in
    suffix (str): anything else you want to add to the log file name

    Returns
    -------
    str
        The file path to the log file
    """
    try:
        # PRECONDITIONS
        Path(folder).mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise OSError(f"Error creating directory '{folder}': {e}")
    # MAIN FUNCTION
    log_filename = get_time(incl_timezone=False) + "_" + suffix + ".log"
    log_filepath = os.path.join(folder, log_filename)

    return log_filepath


def init_log(
    filename: str, display: bool = False, logger_id: str | None = None
) -> logging.Logger:
    """
    - Custom python logger configuration (basicConfig())
        with two handlers (for stdout and for file)
    - from: https://stackoverflow.com/a/44760039
    - Keeps a log record file of the python application, with option to
        display in stdout

    Parameters
    ----------
    filename (str): filepath to log record file
    - display (bool): whether to print the logs to whatever standard output
    - logger_id (str): an optional identifier for yourself,
        if None then defaults to 'root'

    Returns
    -------
    logging.Logger
        The logger object

    Examples
    -----
    >>> logger = init_log('logs/tmp.log', display=True)
    >>> logger.info('Loading things')
    [2023-10-20 10:38:03,074] root: INFO - Loading things
    """
    # PRECONDITIONALS
    assert isinstance(filename, str), "Filename must be a string"
    assert (
        isinstance(logger_id, str) or logger_id is None
    ), "logger_id must be a string or None"

    # MAIN FUNCTION
    # init handlers
    file_handler = logging.FileHandler(filename=filename)
    stdout_handler = logging.StreamHandler(stream=sys.stdout)
    if display:
        handlers = [file_handler, stdout_handler]
    else:
        handlers = [file_handler]

    # instantiate the logger
    logger = logging.getLogger(logger_id)
    logger.setLevel(logging.DEBUG)
    # logger configuration
    # ! logging.basicConfig has no effect if called once anywhere in the code
    # ! set handlers and format for the logger manually
    # Reset any existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Set up the new handlers and format
    formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s - %(message)s")
    for handler in handlers:
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logging.getLogger("matplotlib.font_manager").disabled = True

    return logger


def get_logger(
    log_suffix, folder="logs", display=True, logger_id="metagenomics_parser"
) -> tuple[logging.Logger, str]:
    """
    Initialize the logger with a log file name that includes an optional suffix.

    Parameters
    ----------
    log_suffix : str
        A string to append to the log file name.

    Returns
    -------
    tuple[logging.Logger, str]
        A tuple containing the logger instance and the log file path.
    """
    # Generate log file name
    log_file = generate_log_filename(folder=folder, suffix=log_suffix)

    # Initialize logger
    logger = init_log(log_file, display=display, logger_id=logger_id)

    # Log the path to the log file
    logger.info(f"Path to log file: {log_file}")

    return logger, log_file
