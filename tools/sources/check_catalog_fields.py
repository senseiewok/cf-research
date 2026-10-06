"""Check that every field used in sources/catalog.yaml is documented in the Fields table of sources/README.md.

Exit 0 when every used key is documented (a documented but unused key is only a warning), 1 when a used key is undocumented,
2 on an input problem. Standard library plus PyYAML; no network. Run from anywhere: the default paths come from this file's location.
"""
import argparse
import re
import sys
from pathlib import Path
import yaml


def parse_args():
    parser = argparse.ArgumentParser(description="Check catalog fields are documented in README.")
    parser.add_argument("--catalog", default=None, help="Path to catalog YAML file")
    parser.add_argument("--readme", default=None, help="Path to README markdown file")
    return parser.parse_args()


def get_default_paths():
    root = Path(__file__).resolve().parents[2]
    return root / "sources" / "catalog.yaml", root / "sources" / "README.md"


def extract_used_keys(catalog_path):
    try:
        with open(catalog_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Catalog file not found: {catalog_path}")
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML in catalog: {e}")

    if not isinstance(data, dict):
        raise ValueError("Catalog must be a mapping with 'sources' list")

    sources = data.get("sources")
    if sources is None or not isinstance(sources, list):
        raise ValueError("Catalog has no 'sources' list")

    used_keys = set()
    for item in sources:
        if isinstance(item, dict):
            used_keys.update(item.keys())

    defaults = data.get("defaults")
    if isinstance(defaults, dict):
        used_keys.update(defaults.keys())

    return used_keys


def extract_documented_keys(readme_path):
    try:
        with open(readme_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"README file not found: {readme_path}")

    # Find the Fields table: header row is exactly `| Field | Meaning |` (ignoring spacing around pipes)
    in_fields_table = False
    seen_separator = False
    documented_keys = set()
    past_table = False

    for line in lines:
        stripped = line.strip()

        if not in_fields_table:
            # Check if this is the Fields table header
            # The header must be exactly `| Field | Meaning |` ignoring spacing around pipes
            # So we normalize: remove spaces around pipes
            normalized = re.sub(r'\s*\|\s*', '|', stripped)
            if normalized == '|Field|Meaning|':
                in_fields_table = True
                seen_separator = False
                continue

        if in_fields_table:
            if not seen_separator:
                # Next line should be the separator row
                # Separator rows typically look like |---|---| or | --- | --- |
                if re.match(r'^\s*\|[\s\-:]+\|\s*$', stripped) or (
                    '|' in stripped and '-' in stripped and not any(c.isalpha() for c in stripped.replace('|', '').replace('-', '').replace(':', '').strip())
                ):
                    seen_separator = True
                    continue
                else:
                    # Not a valid separator, stop looking
                    in_fields_table = False
                    past_table = True
                    continue

            # Now we're reading data rows until a line that doesn't start with |
            if not stripped.startswith('|'):
                break  # Stop reading at first line after table that doesn't start with |

            # Extract words in backticks from the FIRST cell
            # Split by pipes, take the second element (first cell is empty due to leading pipe)
            parts = stripped.split('|')
            # After splitting "| `a`, `b` | meaning |", we get ['', ' `a`, `b` ', ' meaning ', '']
            if len(parts) >= 2:
                first_cell = parts[1]
                # Find all backticked words
                matches = re.findall(r'`([^`]+)`', first_cell)
                for m in matches:
                    documented_keys.add(m.strip())

    if not documented_keys and not past_table:
        raise ValueError("README has no Fields table")

    # We need to verify we actually found the table
    # If we never entered the table, it's an error
    return documented_keys


def main():
    args = parse_args()

    default_catalog, default_readme = get_default_paths()
    catalog_path = Path(args.catalog) if args.catalog else default_catalog
    readme_path = Path(args.readme) if args.readme else default_readme

    try:
        used_keys = extract_used_keys(catalog_path)
    except (FileNotFoundError, ValueError) as e:
        print(f"error: {e}")
        return 2

    try:
        documented_keys = extract_documented_keys(readme_path)
    except (FileNotFoundError, ValueError) as e:
        print(f"error: {e}")
        return 2

    undocumented = used_keys - documented_keys
    unused = documented_keys - used_keys

    if undocumented:
        sorted_undoc = sorted(undocumented)
        print("undocumented: " + ", ".join(sorted_undoc))
        # Spec: the unused warning is printed whenever there are unused keys, whether or not anything is undocumented.
        if unused:
            sorted_unused = sorted(unused)
            print("unused (warning): " + ", ".join(sorted_unused))
        return 1

    if unused:
        sorted_unused = sorted(unused)
        print("unused (warning): " + ", ".join(sorted_unused))

    n = len(used_keys)
    print(f"fields ok: {n} keys documented")
    return 0


if __name__ == "__main__":
    sys.exit(main())
