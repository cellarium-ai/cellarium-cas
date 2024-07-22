import platform
import subprocess
import typing as t

from typing_extensions import Self


def find_and_kill_process(port: int, verbose: bool = False):
    """
    Find and kill processes using a specific port.

    :param port: The port number to search for and kill processes.
    :type port: int

    :raises OSError: If the operating system is not supported.

    :return: None
    """

    system = platform.system()

    if system == "Linux" or system == "Darwin":  # macOS is Darwin
        find_command = f"lsof -i :{port}"
        kill_command = "kill -9"

    elif system == "Windows":
        find_command = f"netstat -ano | findstr :{port}"
        kill_command = "taskkill /PID"

    else:
        raise OSError(f"Unsupported operating system: {system}")

    try:
        # Find the process using the port
        result = subprocess.run(find_command, shell=True, check=True, capture_output=True, text=True)
        output = result.stdout

        if system == "Linux" or system == "Darwin":
            lines = output.strip().split("\n")
            for line in lines[1:]:
                pid = int(line.split()[1])
                subprocess.run(f"{kill_command} {pid}", shell=True)
                if verbose:
                    print(f"Killed process with PID: {pid}")

        elif system == "Windows":
            lines = output.strip().split("\n")
            for line in lines:
                pid = int(line.split()[-1])
                subprocess.run(f"{kill_command} {pid} /F", shell=True)
                if verbose:
                    print(f"Killed process with PID: {pid}")

        if verbose:
            print(f"All processes using port {port} have been killed.")

    except subprocess.CalledProcessError as e:
        if verbose:
            print(f"No process found using port {port}. {e}")


T = t.TypeVar("T")


class ConfigValue(t.Generic[T]):
    """
    Class to store a configuration value and its original value.  Useful for tracking uncommitted changes.

    Note: there is no typechecking to ensure that the original value and the new value are of the same type.
    """

    def __init__(self, value: T):
        self.__value__ = value
        self.__original_value__ = value
        self.__dirty_value__ = None
        self.__is_dirty__ = False

    def set(self, value: T) -> Self:
        self.__dirty_value__ = value
        self.__is_dirty__ = True
        return self

    def rollback(self) -> Self:
        self.__dirty_value__ = None
        self.__is_dirty__ = False
        return self

    def commit(self) -> Self:
        if self.__is_dirty__:
            self.__value__ = self.__dirty_value__
        self.__is_dirty__ = False
        return self

    def get(self, dirty_read: bool = False) -> T:
        if dirty_read and self.__is_dirty__:
            return self.__dirty_value__
        return self.__value__

    def is_dirty(self) -> bool:
        return self.__is_dirty__

    def reset(self) -> Self:
        self.__value__ = self.__original_value__
        self.__dirty_value__ = None
        self.__is_dirty__ = False
        return self
