#!/usr/bin/env python3
"""
Extract cross-sections from MadGraph banner files and create a table.

This script:
1. Scans the <OUTPUT_FOLDER>/results/ directory for Events_* directories
2. Parses the directory names to extract grid point parameters
3. Reads the cross-section from the banner.txt files
4. Creates a formatted table output file
"""

import re
import sys
from pathlib import Path

def parse_directory_name(dirname):
    """
    Parse directory name like 'Events_st0p35_tb3_mA200_ma50'
    Returns a dictionary with sintheta, tanbeta, mA, ma
    """
    # Remove 'Events_' prefix
    name = dirname.replace('Events_', '')

    # Extract sintheta (st0p35 -> 0.35)
    st_match = re.search(r'st([0-9]+)p([0-9]+)', name)
    if st_match:
        sintheta = float(f"{st_match.group(1)}.{st_match.group(2)}")
    else:
        sintheta = None

    # Extract tanbeta (tb3 -> 3)
    tb_match = re.search(r'tb([0-9]+)', name)
    if tb_match:
        tanbeta = float(tb_match.group(1))
    else:
        tanbeta = None

    # Extract mA (mA200 -> 200)
    mA_match = re.search(r'mA([0-9]+)', name)
    if mA_match:
        mA = float(mA_match.group(1))
    else:
        mA = None

    # Extract ma (ma50 -> 50)
    ma_match = re.search(r'ma([0-9]+)', name)
    if ma_match:
        ma = float(ma_match.group(1))
    else:
        ma = None

    return {
        'sintheta': sintheta,
        'tanbeta': tanbeta,
        'mA': mA,
        'ma': ma,
        'dirname': dirname
    }

def extract_cross_section(banner_path):
    """
    Extract cross-section from banner file.
    Looks for line: "#  Integrated weight (pb)  :       0.00017319891580000002"
    """
    try:
        with open(banner_path, 'r') as f:
            for line in f:
                if 'Integrated weight (pb)' in line:
                    # Extract the number after the colon
                    match = re.search(r':\s*([\+\-\d.eE]+)', line)
                    if match:
                        return float(match.group(1))
    except FileNotFoundError:
        return None
    except Exception as e:
        print(f"Error reading {banner_path}: {e}")
        return None

    return None

def find_banner_file(results_dir, event_dir):
    """
    Find the banner file in the results directory.
    Looks for: results/Events_*/run_*/run_*_tag_*_banner.txt
    """
    event_path = Path(results_dir) / event_dir

    # Look for run_* directories
    run_dirs = sorted(event_path.glob('run_*'))

    for run_dir in run_dirs:
        # Look for banner files
        banner_files = sorted(run_dir.glob('*_banner.txt'))
        if banner_files:
            return str(banner_files[0])

    return None

def main():
    if len(sys.argv) != 2:
        prog = Path(sys.argv[0]).name
        print(f"Usage: {prog} <OUTPUT_FOLDER>")
        print("Example: python3 extract_cross_sections.py bbdm_2HDMa_type1_case1_scan")
        return 2

    base_dir = Path(sys.argv[1]).expanduser()
    results_dir = base_dir / 'results'

    if not results_dir.exists():
        print(f"Error: {results_dir} directory not found!")
        return 2

    # Collect all data
    data = []

    # Find all Events_* directories
    event_dirs = sorted([d.name for d in results_dir.iterdir()
                        if d.is_dir() and d.name.startswith('Events_')])

    print(f"Found {len(event_dirs)} event directories")
    count = 0
    for event_dir in event_dirs:
        # Parse directory name
        params = parse_directory_name(event_dir)

        # Find banner file
        banner_path = find_banner_file(results_dir, event_dir)

        if banner_path:
            cross_section = extract_cross_section(banner_path)
            if cross_section is not None:
                params['cross_section'] = cross_section
                data.append(params)
                count += 1
                if count % 100 == 0:
                    print(f"Processed {count}th event: {event_dir} -> {cross_section:.15f} pb")
            else:
                print(f"Warning: Could not extract cross-section from {banner_path}")
        else:
            print(f"Warning: Could not find banner file for {event_dir}")

    # Sort data by parameters for better readability
    data.sort(key=lambda x: (
        x.get('sintheta') or 0,
        x.get('tanbeta') or 0,
        x.get('mA') or 0,
        x.get('ma') or 0
    ))

    # Write table to file
    output_file = base_dir / 'cross_section_table.txt'

    with open(output_file, 'w') as f:
        # Write header
        f.write("# Cross-section table\n")
        f.write("# Format: sintheta | tanbeta | mA | ma | cross_section (pb)\n")
        f.write("#" + "="*80 + "\n")
        f.write(f"{'sintheta':<12} {'tanbeta':<10} {'mA':<8} {'ma':<8} {'cross_section (pb)':<20}\n")
        f.write("-" * 80 + "\n")

        # Write data
        for entry in data:
            sintheta = entry.get('sintheta', 'N/A')
            tanbeta = entry.get('tanbeta', 'N/A')
            mA = entry.get('mA', 'N/A')
            ma = entry.get('ma', 'N/A')
            cross_section = entry.get('cross_section', 'N/A')

            if isinstance(cross_section, float):
                f.write(f"{sintheta:<12.4f} {tanbeta:<10.1f} {mA:<8.0f} {ma:<8.0f} {cross_section:<20.15f}\n")
            else:
                f.write(f"{sintheta:<12} {tanbeta:<10} {mA:<8} {ma:<8} {cross_section:<20}\n")

    print(f"\nTable written to: {output_file}")
    print(f"Total entries: {len(data)}")

    return 0

if __name__ == '__main__':
    raise SystemExit(main())
