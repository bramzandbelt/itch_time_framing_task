# !/Users/bramzandbelt/anaconda/envs/psychopyenv/bin python
# -*- coding: utf-8 -*-

"""
This script runs the inertemporal choice task

"""

__author__      = "bramzandbelt (bramzandbelt@gmail.com)"
__copyright__   = "Copyright (c) 2018 Bram Zandbelt"
__license__     = "CC BY 4.0"
__version__     = "0.1"
__vcs_id__      = ""

###############################################################################
# IMPORT MODULES
###############################################################################
import psychopy as pp
from psychopy import iohub as ppiohub
from psychopy import gui as ppgui
import itch as itch
import sys
import glob

# Assure that text files are read in as Unicode (utf8) instead of ASCII
reload(sys)
sys.setdefaultencoding('utf8')


class ExperimentRuntime(ppiohub.ioHubExperimentRuntime):
    """

    Create class

    """

    def run(self,*args):
        """

        The run method contains your experiment logic. It is equal to what would be in your main psychopy experiment
        script.py file in a standard psychopy experiment setup. That is all there is too it really.

        :param args:

        """
        # from psychopy import gui

        # 1. Parse and process configuration
        # ---------------------------------------------------------------------
        mod_dir = args[0][1]
        config  = itch.init_config(self,mod_dir)
        session_ix = config['session']['session_ix']

        # 2. Run practice blocks
        # ---------------------------------------------------------------------
        # if config['practice']['enable']:
        #
        #     p_trial_list = config['practice']['trial_list']
        #     p_trial_list = p_trial_list[p_trial_list.session_ix == session_ix]
        #
        #     itch.run_stage(config = config,
        #               stage_id = 'practice',
        #               trial_list = p_trial_list)

        # 2. Run indifference point procedure
        # ---------------------------------------------------------------------
        # itch.run_ip_procedure(config)
        # pp.core.quit()


        # 5. Run experimental blocks
        # ---------------------------------------------------------------------
        if config['experiment']['enable']:

            e_trial_list = config['experiment']['trial_list']
            e_trial_list = e_trial_list[e_trial_list.session_ix == session_ix]

            itch.run_stage(config = config,
                           stage_id = 'experiment',
                           trial_list = e_trial_list)

        # 6.Terminate experiment
        # ---------------------------------------------------------------------
        itch.present_instruction(config, 'end')
        pp.core.quit()

####### Main Script Launching Code Below #######

if __name__ == "__main__":
    import os
    import psychopy as pp

    def main(mod_dir):
        """

        :param mod_dir: Directory where module files reside

        """

        config_dir = os.path.normcase(os.path.join(mod_dir,'config/'))

        # Let user select response device
        # ---------------------------------------------------------------------
        rd_config_files = {'Keyboard':
                               'iohub_keyboard.yaml',
                           'Serial':
                               'iohub_serial.yaml',
                           'fORP':
                               'iohub_forp.yaml'}

        info = {'Response Device = ': ['Select',
                                       'Keyboard',
                                       'Serial',
                                       'fORP']}

        # dlg_info = dict(info)
        # info_dlg = ppgui.DlgFromDict(dictionary = dlg_info,
        #                               title = 'Select response device')
        # if not info_dlg.OK:
        #     return -1
        #
        # while dlg_info.values()[0] == u'Select' and info_dlg.OK:
        #         dlg_info = dict(info)
        #         info_dlg = ppgui.DlgFromDict(dictionary=dlg_info,
        #                                       title='SELECT Response device to continue...')
        #
        # if not info_dlg.OK:
        #     return -1

        # Determine which iohub base configuration file to use
        # ---------------------------------------------------------------------
        # os.chdir(config_dir)
        # iohub_base_files = ['Select']
        # iohub_base_files.extend(glob.glob('./iohub_base*.yaml'))
        # os.chdir(mod_dir)
        #
        # iohub_base_config_info = {'IOHub base config file = ': iohub_base_files}
        # iohub_base_config_info = dict(iohub_base_config_info)
        #
        # iohub_base_config_dlg = ppgui.DlgFromDict(dictionary = iohub_base_config_info,
        #                                            title = 'Select IOHub base config file')
        # if not iohub_base_config_dlg.OK:
        #     return -1
        #
        # while iohub_base_config_info.values()[0] == u'Select' and iohub_base_config_dlg.OK:
        #     iohub_base_config_info = dict(iohub_base_config_info)
        #     iohub_base_config_dlg = ppgui.DlgFromDict(dictionary = iohub_base_config_info,
        #                                                title = 'SELECT IOHub base config file to continue...')

        # Merge iohub configuration files
        # ---------------------------------------------------------------------
        base_config_file = os.path.normcase(
            os.path.join(
                config_dir,
                'iohub_base_itch_experiment_macbook.yaml'.replace('./', '')
                # iohub_base_config_info.values()[0].replace('./','')
            )
        )

        resp_dev_config_file = \
            os.path.normcase(
                os.path.join(
                    config_dir,
                    'iohub_keyboard.yaml'
                    # rd_config_files[dlg_info.values()[0]]
                )
            )

        combined_config_file = \
            os.path.normcase(os.path.join(config_dir,
                                          'iohub_config.yaml'
                                          )
                             )

        ExperimentRuntime.mergeConfigurationFiles(base_config_file,
                                                  resp_dev_config_file,
                                                  combined_config_file)

        # Determine which experiment configuration file to use
        # ---------------------------------------------------------------------
        # os.chdir(config_dir)
        # expt_files = ['Select']
        # expt_files.extend(glob.glob('./expt_*.yaml'))
        # os.chdir(mod_dir)
        #
        # expt_config_info = {'Experiment config file = ': expt_files}
        # expt_config_info = dict(expt_config_info)
        #
        # expt_config_dlg = ppgui.DlgFromDict(dictionary = expt_config_info,
        #                                      title = 'Select experiment config file')
        # if not expt_config_dlg.OK:
        #     return -1
        #
        # while expt_config_info.values()[0] == u'Select' and expt_config_dlg.OK:
        #     expt_config_info = \
        #         dict(expt_config_info)
        #     expt_config_dlg = \
        #         ppgui.DlgFromDict(dictionary=expt_config_info,
        #                            title='SELECT experiment config file to continue...')
        #
        # if not expt_config_dlg.OK:
        #     return -1

        # Start the experiment
        # ---------------------------------------------------------------------
        runtime = ExperimentRuntime(config_dir, 'expt_itch.yaml')
        runtime.start(('iohub_keyboard.yaml', mod_dir))
        # runtime = ExperimentRuntime(config_dir, expt_config_info.values()[0])
        # runtime.start((dlg_info.values()[0],mod_dir))

    # Get the current directory, using a method that does not rely on __FILE__
    # or the accuracy of the value of __FILE__.
    mod_dir = ppiohub.module_directory(main)

    # Run the main function, which starts the experiment runtime
    main(mod_dir)