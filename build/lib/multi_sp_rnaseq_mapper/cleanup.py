
import os

def delete_fastq(output_dir):
    """Delete fastq files after alignment."""

    fastq_list = [ f for f in os.listdir(output_dir) if f.endswith(".fastq") ]
    for f in fastq_list:
        os.remove(os.path.join(output_dir, f))

    # try:
    #     if layout == 'paired':
    #         print(f"Deleting paired-end fastq files for {sra_id}...")
    #         fastq1 = os.path.join(output_dir, f"{sra_id}_1.fastq")
    #         fastq2 = os.path.join(output_dir, f"{sra_id}_2.fastq")
    #         os.remove(fastq1)
    #         os.remove(fastq2)
    #     elif layout == 'single':
    #         print(f"Deleting single-end fastq file for {sra_id}...")
    #         single_fastq = os.path.join(output_dir, f"{sra_id}.fastq")
    #         os.remove(single_fastq)
    # except FileNotFoundError as e:
    #     print(f"Error deleting fastq files: {e}")
