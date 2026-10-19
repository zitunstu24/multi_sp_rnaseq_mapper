import os
import subprocess

def download_sra(sra_id, output_dir, layout):
    """Download SRA files using fasterq-dump based on layout (single or paired)."""
    try:
        print(f"Downloading SRA sample {sra_id} as {layout}-end...")
        if layout == 'paired':
            subprocess.run(['fasterq-dump', '--split-files', '-O', output_dir, sra_id], check=True)
        elif layout == 'single':
            subprocess.run(['fasterq-dump', '-O', output_dir, sra_id], check=True)
        else:
            raise ValueError(f"Unknown layout type: {layout}. Should be 'single' or 'paired'.")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading {sra_id}: {e}")
