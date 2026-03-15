"""Profile system for specimux parameter presets.

Profiles allow users to save and reuse parameter configurations for different
workflows (e.g., different primer sets, length filtering, threading).

Profile resolution order:
1. User profiles in ~/.config/specimux/profiles/
2. Bundled profiles in package

Override order: defaults -> profile -> explicit CLI arguments
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Any, Optional, List, Set
import logging
import os
import re
import sys

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore

try:
    from specimux import __version__
except ImportError:
    __version__ = "dev"

# Use importlib.resources for Python 3.9+, fall back to pkg_resources
HAS_IMPORTLIB_RESOURCES = False
HAS_PKG_RESOURCES = False

try:
    from importlib.resources import files as importlib_files
    from importlib.resources import as_file
    HAS_IMPORTLIB_RESOURCES = True
except ImportError:
    try:
        import pkg_resources
        HAS_PKG_RESOURCES = True
    except ImportError:
        pass

logger = logging.getLogger(__name__)

# XDG-compliant config path
XDG_CONFIG_HOME = Path.home() / ".config"
PROFILES_DIR = XDG_CONFIG_HOME / "specimux" / "profiles"

# Valid keys for the specimux tool section (drawn from CLI options)
VALID_SPECIMUX_KEYS = {
    "trim",
    "dereplicate",
    "search-len",
    "index-edit-distance",
    "primer-edit-distance",
    "min-length",
    "max-length",
    "threads",
    "disable-prefilter",
    "disable-preorient",
    "sample-topq",
    "diagnostics",
}


class ProfileError(Exception):
    """Error loading or applying a profile."""
    pass


class ProfileVersionError(ProfileError):
    """Profile version is incompatible with current specimux version."""
    pass


class ProfileValidationError(ProfileError):
    """Profile contains invalid keys."""
    pass


@dataclass
class Profile:
    """A parameter profile for specimux."""
    name: str
    version: str  # e.g., "0.7.*"
    description: str
    specimux: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def load(cls, name: str, check_version: bool = True) -> 'Profile':
        """Load profile by name from user dir or bundled.

        Args:
            name: Profile name (without .yaml extension)
            check_version: If True, validate version compatibility

        Returns:
            Loaded Profile instance

        Raises:
            ProfileError: If profile not found
            ProfileVersionError: If version incompatible
            ProfileValidationError: If profile contains invalid keys
        """
        if yaml is None:
            raise ProfileError(
                "PyYAML is required for profile support. "
                "Install with: pip install pyyaml"
            )

        # Initialize user profiles directory with example on first use
        ensure_user_profiles_dir()

        # Try user profile first
        user_path = PROFILES_DIR / f"{name}.yaml"
        if user_path.exists():
            return cls._load_from_path(user_path, name, check_version)

        # Fall back to bundled profile
        bundled_path = get_bundled_profile_path(name)
        if bundled_path is not None:
            return cls._load_from_path(bundled_path, name, check_version)

        # Profile not found - provide helpful error
        available = list_profiles()
        if available:
            raise ProfileError(
                f"Profile '{name}' not found. Available profiles: {', '.join(available)}"
            )
        else:
            raise ProfileError(
                f"Profile '{name}' not found and no profiles are available. "
                f"Check that profiles are installed in {PROFILES_DIR}"
            )

    @classmethod
    def _load_from_path(cls, path: Path, name: str, check_version: bool) -> 'Profile':
        """Load profile from a specific path."""
        try:
            with open(path, 'r') as f:
                data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ProfileError(f"Invalid YAML in profile '{name}': {e}")
        except IOError as e:
            raise ProfileError(f"Cannot read profile '{name}': {e}")

        if not isinstance(data, dict):
            raise ProfileError(f"Profile '{name}' must be a YAML mapping")

        # Extract fields
        version = data.get('specimux-version', '*')
        description = data.get('description', '')
        specimux = data.get('specimux', {}) or {}

        # Validate version compatibility
        if check_version and not check_version_compatible(version, __version__):
            raise ProfileVersionError(
                f"Profile '{name}' requires specimux version {version}, "
                f"but you have {__version__}.\n\n"
                f"This profile may use parameters that have changed or been removed.\n"
                f"Please update the profile for your version, or copy the bundled\n"
                f"'{name}' profile which is compatible with your version:\n\n"
                f"  cp {PROFILES_DIR}/{name}.yaml {PROFILES_DIR}/{name}.yaml.bak\n"
                f"  specimux --list-profiles  # Will show available profiles"
            )

        # Validate keys (strict validation)
        _validate_profile_keys(name, specimux)

        return cls(
            name=name,
            version=version,
            description=description,
            specimux=specimux,
        )


def check_version_compatible(profile_version: str, current_version: str) -> bool:
    """Check if profile version pattern matches current version.

    Supports wildcards:
    - "0.7.*" matches "0.7.0", "0.7.1", etc.
    - "0.7.0" matches only "0.7.0"
    - "0.*" matches any 0.x release
    - "*" matches any version

    Args:
        profile_version: Version pattern from profile (e.g., "0.7.*")
        current_version: Current specimux version (e.g., "0.7.2")

    Returns:
        True if versions are compatible
    """
    if profile_version == '*':
        return True

    # Convert wildcard pattern to regex
    # Escape dots and convert * to .*
    pattern = profile_version.replace('.', r'\.').replace('*', r'.*')
    pattern = f'^{pattern}$'

    try:
        return bool(re.match(pattern, current_version))
    except re.error:
        # Invalid regex, treat as literal match
        return profile_version == current_version


def _validate_profile_keys(
    name: str,
    specimux: Dict[str, Any],
) -> None:
    """Validate that profile only contains known keys.

    Raises ProfileValidationError for unknown keys.
    """
    unknown = set(specimux.keys()) - VALID_SPECIMUX_KEYS
    if unknown:
        raise ProfileValidationError(
            f"Profile '{name}' contains unknown keys:\n"
            f"  specimux: {', '.join(sorted(unknown))}\n\n"
            f"This may indicate a typo or an option that has been removed.\n"
            f"Please check the profile and fix or remove the invalid keys."
        )


def get_bundled_profile_path(name: str) -> Optional[Path]:
    """Get path to bundled profile using importlib.resources.

    Args:
        name: Profile name (without .yaml extension)

    Returns:
        Path to bundled profile file, or None if not found
    """
    if HAS_IMPORTLIB_RESOURCES:
        try:
            # Python 3.9+ style
            profiles_pkg = importlib_files('specimux.profiles')
            profile_file = profiles_pkg.joinpath(f'{name}.yaml')
            # Check if file exists using as_file context manager
            with as_file(profile_file) as path:
                if path.exists():
                    return path
        except (TypeError, FileNotFoundError, ModuleNotFoundError):
            pass

    if HAS_PKG_RESOURCES:
        try:
            # Fall back to pkg_resources for older Python
            resource_path = pkg_resources.resource_filename(
                'specimux.profiles', f'{name}.yaml'
            )
            path = Path(resource_path)
            if path.exists():
                return path
        except (FileNotFoundError, ModuleNotFoundError):
            pass

    # Last resort: check relative to this file (profiles/__init__.py)
    bundled_dir = Path(__file__).parent
    bundled_path = bundled_dir / f'{name}.yaml'
    if bundled_path.exists():
        return bundled_path

    return None


def list_bundled_profiles() -> List[str]:
    """List names of bundled profiles."""
    profiles = []

    # Try importlib.resources first
    if HAS_IMPORTLIB_RESOURCES:
        try:
            profiles_pkg = importlib_files('specimux.profiles')
            for item in profiles_pkg.iterdir():
                if str(item).endswith('.yaml'):
                    name = Path(str(item)).stem
                    profiles.append(name)
            if profiles:
                return sorted(profiles)
        except (TypeError, FileNotFoundError, ModuleNotFoundError):
            pass

    # Fall back to checking directory relative to this file
    bundled_dir = Path(__file__).parent
    if bundled_dir.exists():
        for yaml_file in bundled_dir.glob('*.yaml'):
            profiles.append(yaml_file.stem)

    return sorted(profiles)


def list_profiles() -> List[str]:
    """List available profiles (user + bundled).

    Returns list of profile names (without .yaml extension).
    User profiles take precedence over bundled profiles with same name.
    """
    profiles: Set[str] = set()

    # User profiles
    if PROFILES_DIR.exists():
        for yaml_file in PROFILES_DIR.glob('*.yaml'):
            profiles.add(yaml_file.stem)

    # Bundled profiles
    profiles.update(list_bundled_profiles())

    return sorted(profiles)


def ensure_user_profiles_dir() -> Path:
    """Ensure user profiles directory exists with example profile.

    On first use, creates the directory and copies an example profile
    to help users create their own profiles.

    This function is safe to call from parallel processes - it uses
    atomic file operations to avoid race conditions.

    Returns:
        Path to user profiles directory
    """
    import tempfile

    PROFILES_DIR.mkdir(parents=True, exist_ok=True)

    example_path = PROFILES_DIR / "example.yaml"

    # Skip if example already exists (common case, avoid extra work)
    if example_path.exists():
        return PROFILES_DIR

    # Skip if user already has other profiles (they don't need the example)
    if any(PROFILES_DIR.glob('*.yaml')):
        return PROFILES_DIR

    # Copy example profile atomically (safe for parallel invocations)
    bundled_example = get_bundled_profile_path('example')
    if bundled_example is not None:
        try:
            # Write to temp file in same directory, then atomic rename
            fd, temp_path = tempfile.mkstemp(
                dir=PROFILES_DIR,
                prefix='.example.',
                suffix='.yaml.tmp'
            )
            try:
                with open(bundled_example, 'rb') as src:
                    os.write(fd, src.read())
            finally:
                os.close(fd)

            # Atomic rename - if file exists, this either succeeds or fails cleanly
            try:
                os.rename(temp_path, example_path)
                logger.info(f"Created example profile at {example_path}")
            except OSError:
                # Another process won the race - that's fine
                if os.path.exists(temp_path):
                    os.unlink(temp_path)
        except Exception as e:
            # Non-fatal - profile system works without example in user dir
            logger.debug(f"Could not create example profile: {e}")

    return PROFILES_DIR


def apply_profile_to_args(
    args,
    profile: Profile,
    tool: str,
    explicit_args: Set[str]
) -> None:
    """Apply profile values to args, respecting explicit CLI overrides.

    Args:
        args: argparse Namespace to modify
        profile: Profile to apply
        tool: Must be 'specimux'
        explicit_args: Set of argument names that were explicitly provided on CLI
    """
    if tool == 'specimux':
        profile_values = profile.specimux
    else:
        raise ValueError(f"Unknown tool: {tool}")

    for key, value in profile_values.items():
        # Convert YAML key (with dashes) to argparse attribute name (with underscores)
        attr_name = key.replace('-', '_')

        # Only apply if not explicitly set on command line
        if attr_name not in explicit_args:
            if hasattr(args, attr_name):
                logger.debug(f"Profile '{profile.name}': setting {attr_name}={value}")
                setattr(args, attr_name, value)
            else:
                # This shouldn't happen if validation passed, but log it
                logger.warning(
                    f"Profile '{profile.name}': unknown attribute '{attr_name}'"
                )


def print_profiles_list(tool: str = 'specimux') -> None:
    """Print available profiles to stdout.

    Args:
        tool: Tool name for usage hint (default: 'specimux')
    """
    if yaml is None:
        print("Profile support requires PyYAML. Install with: pip install pyyaml")
        return

    # Initialize user profiles directory with example on first use
    ensure_user_profiles_dir()

    profiles = list_profiles()

    if not profiles:
        print(f"No profiles found.")
        print(f"\nProfiles are stored in: {PROFILES_DIR}")
        return

    print(f"Available profiles:\n")

    for name in profiles:
        try:
            # Load without version check to show all profiles
            profile = Profile.load(name, check_version=False)

            # Check if it's a user profile or bundled
            user_path = PROFILES_DIR / f"{name}.yaml"
            source = "user" if user_path.exists() else "bundled"

            # Check version compatibility
            compatible = check_version_compatible(profile.version, __version__)
            compat_str = "" if compatible else " [INCOMPATIBLE]"

            print(f"  {name} ({source}){compat_str}")
            if profile.description:
                print(f"    {profile.description}")
            print(f"    Version: {profile.version}")
            print()

        except ProfileError as e:
            print(f"  {name} [ERROR: {e}]")
            print()

    print(f"Usage: {tool} -p <profile> [other options]")
    print(f"Profile directory: {PROFILES_DIR}")
