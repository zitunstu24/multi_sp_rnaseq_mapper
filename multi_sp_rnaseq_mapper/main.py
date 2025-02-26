import argparse
from multi_sp_rnaseq_mapper.process_samples import process_samples
from multi_sp_rnaseq_mapper.DEGs import get_DEGs

def main():
    parser = argparse.ArgumentParser(description="SRA Alignment Tool")
    parser.add_argument('--config', required=True, help="Path to the config YAML file")

    args = parser.parse_args()

    # Process samples based on config
    process_samples(args.config)
    get_DEGs(args.config)


if __name__ == "__main__":
    main()
