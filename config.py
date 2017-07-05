import ConfigParser

class MCSNConfig:
    """
    A container class of the configuration data.
    Saves the configuration data as a Python dictionary.
    """
    def __init__(self, file_path=None):
        # The data dict has option values as Python objects,
        # some of them being strings.
        self.data = {}
        self.parser = ConfigParser.SafeConfigParser()
        # attribute options has the following structure:
        # options['section']['option'] = 'label'
        self.options = {
            'Network Data': {
                'streets': 'Street Names',
                'branch_points': 'Branch Points',
                'joints': 'Joints',
                'homology_classes': 'Homology Classes',
            },
            'Computation Parameters': {
                'iterations': 'Number of iterations',
            },
        }

        # # {deprecated option: option,}
        # self.deprecated_options = {
        #     'phase_range': 'phase',
        #     'punctures': 'irregular_punctures',
        #     'size_of_neighborhood': 'size_of_bp_neighborhood',
        #     'size_of_bin': None,
        #     'size_of_ramification_pt_cutoff': None,
        #     'n_processes': None,
        # }

        if file_path is not None:
            print ('Loading configuration from {}...'.format(file_path))
            with open(file_path, 'r') as fp:
                self.read(fp)
            print ('Finished loading configuration from {}.\n\n'
                        .format(file_path))

    def __setitem__(self, option, value):
        try:
            # Update the data dict.
            self.data[option] = value
            # Update the parser
            for section in self.options:
                if option in self.options[section]:
                    self.parser.set(section, option, str(value))
        except KeyError:
            # Not an available the option, raise an error.
            raise KeyError('Unknown option \'{}\'.'.format(option))

    def __getitem__(self, option):
        try:
            return self.data[option]
        except KeyError:
            # Not an available the option, raise an error.
            raise KeyError('Unknown option \'{}\'.'.format(option))

    def __delitem__(self, option):
        try:
            # Remove the option from the data dict.
            del self.data[option]
            # Remove the option from the parser.
            for section in self.options:
                if option in self.options[section]:
                    self.parser.remove_option(section, option)
        except KeyError:
            raise KeyError('Unknown option \'{}\'.'.format(option))

    def get_label(self, option):
        for section in self.options.keys():
            if option in self.options[section]:
                return self.options[section][option]

        raise ValueError('Unknown option \'{}\'.'.format(option))

    def keys(self):
        return self.data.keys()

    def iteritems(self):
        return self.data.iteritems()

    def read(self, config_file):
        """
        Read an .ini file and load the configuration data.
        """
        self.parser.readfp(config_file)

        for section in self.parser.sections():
            try:
                not_configured_options = self.options[section].keys()
            except KeyError:
                self.parser.remove_section(section)
                raise ValueError('Unknown section \'{}\'. in the config file.'
                                 .format(section))

            for option in self.parser.options(section):
                parser_value = self.parser.get(section, option)

                # # Check deprecated options
                # if option in self.deprecated_options:
                #     old_option = option
                #     option = self.deprecated_options[old_option]
                #     logger.warning(
                #         'Option \'{}\' is deprecated.'.format(old_option)
                #     )
                #     if option is None:
                #         continue
                #     logger.warning('Use \'{}\' instead.'.format(option))
                #     self.parser.remove_option(section, old_option)
                #     self.parser.set(section, option, parser_value)

                if (parser_value == 'None'):
                    value = eval(parser_value)
                elif (
                    section == 'Network Data' or 
                    section == 'Computation Parameters'
                ):
                    value = eval(parser_value)
                # elif (section == 'numerical parameters'):
                #     value = eval(parser_value)
                else:
                    raise ValueError(
                        'Option \'{}\' in an unknown section \'{}\'.'
                        .format(option, section)
                    )

                try:
                    not_configured_options.remove(option)
                except ValueError:
                    raise ValueError('Unknown option \'{}\'.'.format(option))

                self.data[option] = value

            for option in not_configured_options:
                self.parser.set(section, option, 'None')
                self[option] = None

    def save(self, file_path=None, logger_name=None,):
        if file_path is None:
            return None
        print('Saving configuration to {}.'.format(file_path))
        with open(file_path, 'w') as fp:
            self.parser.write(fp)
        logger.info('Finished saving configuration to {}.'.format(file_path))

