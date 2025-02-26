
import os
import subprocess

def convert_gff3_to_gtf(gff3_dir, output_dir=None):
    """Convert all .gff3 files in a directory to .gtf format using gffread."""
    
    # Ensure output directory exists
    if output_dir is None:
        output_dir = gff3_dir  # Save to the same directory if no output_dir specified
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop through each .gff3 file in the directory
    for file_name in os.listdir(gff3_dir):
        if file_name.endswith(".gff3"):
            gff3_file = os.path.join(gff3_dir, file_name)
            gtf_file = os.path.join(output_dir, file_name.replace(".gff3", ".gtf"))

            # Convert gff3 to gtf using gffread
            try:
                print(f"Converting {gff3_file} to {gtf_file}...")
                subprocess.run(['gffread', gff3_file, '-T', '-o', gtf_file], check=True)
                print(f"Conversion complete: {gtf_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error converting {gff3_file}: {e}")
