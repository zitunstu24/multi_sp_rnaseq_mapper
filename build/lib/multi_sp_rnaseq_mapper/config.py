import yaml

class Config:
    def __init__(self, config_file='config/default_config.yml'):
        with open(config_file, 'r') as file:
            self.config = yaml.safe_load(file)

    def get(self, section, key, default=None):
        """Get a value from the configuration file."""
        return self.config.get(section, {}).get(key, default)

    def get_master_table(self):
        """Get the path to the master table from the configuration file."""
        return self.config.get('master_table')