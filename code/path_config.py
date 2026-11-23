"""Centralized path loader reading config_paths.yaml from the repository root."""
from pathlib import Path
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = REPO_ROOT / "config_paths.yaml"


def load_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def resolve_path(key: str, subkey: str | None = None) -> Path:
    cfg = load_config()
    val = cfg["paths"][key]
    if subkey is not None:
        val = val[subkey]
    p = Path(val)
    return p if p.is_absolute() else (REPO_ROOT / p).resolve()
