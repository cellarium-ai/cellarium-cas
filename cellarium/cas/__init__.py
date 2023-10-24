from cellarium.cas.client import CASClient  # noqa
import toml

def get_current_version_from_toml(pyproject_path="pyproject.toml"):
    with open(pyproject_path, "r") as file:
        data = toml.load(file)
        version = data.get("tool", {}).get("setuptools", {}).get("version", "")
    return version

print(get_current_version_from_toml(), "ASDASDASDAS")