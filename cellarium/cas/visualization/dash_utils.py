import subprocess
import platform


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
