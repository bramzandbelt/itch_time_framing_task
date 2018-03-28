# !/Users/bramzandbelt/anaconda/envs/psychopyenv/bin python
# -*- coding: utf-8 -*-
__author__ = "bramzandbelt"

# Import packages
# from psychopy import visual, monitors, core, event, iohub, info, gui
import psychopy as pp
from psychopy import info as ppinfo
import numpy as np
import pandas as pd
import calendar
import monthdelta
import os
import random  # for setting random number generator seed
import re
import serial
import time
import itertools as it

from pprint import pprint
from psychopy.hardware.emulator import launchScan

print 'psychopy version: %s' % ppinfo.psychopyVersion
print 'numpy version: %s' % np.__version__
print 'pandas version: %s' % pd.__version__

# This is to enable printing of all data frame columns
pd.set_option('display.max_columns', None)


def collect_response(rd, kb, *args, **kwargs):
    """
    Collect responses from response device and keyboard

    This function is called at two stages in the experiment:
    1. Instruction
        Returns count of response keys, escape keys, and other keys to monitor
    2. Experiment
        Updates the log variable with key count and key time information

    Parameters
    ----------
    rd          : dict
                specifies response device properties

    kb          : dict
                specifies keyboard properties

    log         : pandas.core.frame.DataFrame
                trial log

    other_keys  : list (optional)
                specifies which other keys (from response device or keyboard)
                to monitor

    Returns
    -------
    log         : pandas.core.frame.DataFrame (optional)
                trial log; collect_responses fills in values for key_count and
                key_time variables.

    Usage
    -----
    # For collecting experimental data
    log = collect_response(rd,kb,log)

    # For instruction screens
    key_count = collect_response(rd,kb)

    """
    triggered = None
    other_keys_pressed = None
    other_keys = None
    log = None

    if len(args) == 0:
        if kwargs:
            if 'log' in kwargs:
                log = kwargs.get('log')
            if 'other_keys' in kwargs:
                other_keys = kwargs.get('other_keys')
    elif len(args) == 1:
        other_keys = args[0]
    elif len(args) == 2:
        other_keys = args[0]
        log = args[1]
    elif len(args) > 2:
        # TODO: Add error message
        pass

    # Process inputs
    # -------------------------------------------------------------------------
    rsp_keys = []
    for sublist in rd['settings']['rsp_keys']:
        for item in sublist:
            rsp_keys.append(item)
    key_key = rd['settings']['key_key']
    time_key = rd['settings']['time_key']
    rd_class = rd['settings']['class']
    esc_keys = kb['settings']['esc_keys']

    # Define dynamic variables
    # -------------------------------------------------------------------------
    key_count = {key: 0 for key in rsp_keys}
    key_time = {key: [] for key in rsp_keys}

    # Determine response identity and response
    # -------------------------------------------------------------------------
    rd_events = rd['client'].getEvents()

    for ev in rd_events:

        if rd_class == 'Keyboard':

            # Have any abort keys been pressed?
            if any([re.findall(key, ev.key) for key in esc_keys]):
                print('Warning: Escape key pressed - Experiment is terminated')
                pp.core.quit()

            # Have any other keys been pressed (e.g. toggle keys for moving
            # between instruction screens)
            if other_keys:
                if any([re.findall(key, ev.key) for key in
                        other_keys]):
                    other_keys_pressed = ev.key
                    print other_keys_pressed
                    return key_count, other_keys_pressed

            # If any of the event keys are in event data
            if any([re.findall(key, ev.key) for key in rsp_keys]):
                key_count[ev.key] += 1
                if isinstance(log, pd.DataFrame):
                    # Only log time of first response key event
                    if not any([key_time[key] for key in rsp_keys]):
                        key_time[ev.key] = ev.time - log.iloc[0]['trial_ons']

        elif rd_class == 'Serial':

            for kbev in kb['client'].getEvents():

                # Have any abort keys been pressed?
                if any([re.findall(key, kbev.key) for key in esc_keys]):
                    print(
                        'Warning: Escape key pressed - Experiment is '
                        'terminated')
                    pp.core.quit()

                if other_keys:
                    if any([re.findall(key, kbev.key) for key in
                            other_keys]):
                        other_keys_pressed = kbev.key
                        return key_count, other_keys_pressed

                # If any of the event keys are in event data
                if any([re.findall(key, ev.data) for key in rsp_keys]):
                    key_count[ev.data] += 1

                    # Only log time of first response key event
                    if not any([key_time[key] for key in rsp_keys]):
                        key_time[ev.key] = ev.time - log.iloc[0]['trial_ons']

    # For each key, only response times of first two events are stored
    if isinstance(log, pd.DataFrame):


        # Determine choice
        choice_key, choice_time = \
            next((k, v) for k, v in key_time.items() if v)

        if log.iloc[0]['ll_side'] == 'left':
            key_mapping = dict(zip(rsp_keys, ('ll', 'ss')))
        elif log.iloc[0]['ll_side'] == 'right':
            key_mapping = dict(zip(rsp_keys, ('ss', 'll')))

        # A choice can be made only once all relevant stimuli are onscreen,
        # so determine onset of last relevant stimulus
        last_stim_ons = \
            max([log.iloc[0][col]
                 for col in log.iloc[0].keys()
                 if col in ['choice_instr_ons', 'ss_opt_ons', 'll_opt_ons']])

        # Log events
        for key in rsp_keys:
            log.iloc[0]['key_count_' + key] = key_count[key]
        log.iloc[0]['choice'] = key_mapping[choice_key]
        log.iloc[0]['rt'] = choice_time - last_stim_ons

        return log



    else:

        return key_count, other_keys_pressed

def define_stimulus(window, stim_info, *args):
    # type: (object, object, object) -> object
    """
    Make PsychoPy stimuli based on user input

    Currently, only stimuli of type TextStim and ImageStim are supported.

    Parameters
    ----------
    window      : psychopy.visual.window.Window
                PsychoPy window object, in which stimuli are presented

    stim_info    : dict
                stimulus information specified by user in experiment
                configuration file

    stim_type    : str or unicode (optional)
                stimulus type

    Returns
    -------
    stimulus    : dict
                specifies properties of the stimuli

    """

    # Process inputs
    # -------------------------------------------------------------------------

    # TODO: See if assertions require __debug__ variable
    assert type(window) is pp.visual.window.Window, \
        'window is not of class visual.window.Window'
    assert type(stim_info) is dict, \
        'stim_info is not of type dict'
    for key in ['type', 'name']:
        assert stim_info.has_key(key), \
            'stim_info does not contain key {0:s}'.format(key)

    if len(args) == 1:
        stim_type = args[0]
    else:
        stim_type = ''

    n_stimulus = 1

    # Initialize list of stimuli
    stimulus = [None] * n_stimulus

    for i in range(n_stimulus):

        stimulus[i] = init_stimulus(window,stim_info['type'])

        assert (type(stimulus[i]) is
                pp.visual.text.TextStim or
                pp.visual.image.ImageStim), \
                "stimulus is neither a TextStim nor an ImageStim"

        # Set stimulus name
        stimulus[i].name = stim_type

        # if isinstance(stim_info['name'],list):
        #     stimulus[i].name = ''.join([stim_type,'_',stim_info['name'][i]])
        # else:
        #     stimulus[i].name = ''.join([stim_type,'_',stim_info['name']])

        # Text stimuli
        # ---------------------------------------------------------------------
        if type(stimulus[i]) is pp.visual.text.TextStim:

            # Set stimulus content
            # TODO: Include function here that defines stimulus content
            # based on input


            # Set other parameters
            if 'font_file' in stim_info:
                if stim_info['font_file']:
                    stimulus[i].fontFiles = stim_info['font_file']
            if 'font' in stim_info:
                if stim_info['font']:
                    stimulus[i].setFont(stim_info['font'])
            if 'ori' in stim_info:
                if stim_info['ori']:
                    stimulus[i].setOri(stim_info['ori'])
            if 'height' in stim_info:
                if stim_info['height']:
                    stimulus[i].setHeight(stim_info['height'])
            if 'pos' in stim_info:
                if stim_info['pos']:
                    stimulus[i].setPos(stim_info['pos'])
            if 'color' in stim_info:
                if stim_info['color']:
                    stimulus[i].setColor(stim_info['color'], 'rgb255')
            if 'opacity' in stim_info:
                if stim_info['opacity']:
                    stimulus[i].setOpacity(stim_info['opacity'])

        # Image stimuli
        # ---------------------------------------------------------------------
        else:
            # TODO: Implement UserWarning exception here.
            None


    return stimulus
def evaluate_trial(eval_data,feedback_dur,window,stimuli,log):
    """
    Evaluate trial performance and present trial feedback

    Parameters
    ----------
    eval_data    : dict
                specifies information about how trial should be evaluated,
                this includes:

                eval_data_file    : dict
                                specifies trial evaluation data files for each
                                of the possible response devices
                eval_data        : pandas.core.frame.DataFrame
                                all stimulus-response combinations specified
                                in the trial evaluation data file

                correct         : pandas.core.series.Series
                                accuracy for each of the stimulus-response
                                combinations

                responseType    : pandas.core.series.Series
                                response type for each of the stimulus-response
                                combinations

                feedback        : pandas.core.series.Series
                                feedback for each of the stimulus-response
                                combinations

                trialType       : pandas.core.series.Series
                                trial type for each of the stimulus-response
                                combinations

    feedback_dur : float
                trial feedback duration (in seconds)

    window      : psychopy.visual.window.Window
                PsychoPy window object, in which stimuli are presented

    stimuli     : dict
                specifies PsychoPy stimuli, including the feedback and inter-
                trial interval stimulus

    log         : pandas.core.frame.DataFrame
                trial log

    Returns
    -------
    log         : pandas.core.frame.DataFrame
                trial log; evaluate_trial fills in values for the following
                variables: trialCorrect, trialType, responseType, and
                trialFeedback

    """

    # Process inputs
    # =========================================================================

    # Assertions
    # -------------------------------------------------------------------------

    # specific to eval_data and log (window and stimulis should have been
    # checked already and not have been updated

    # Define dynamic variables
    # -------------------------------------------------------------------------
    # trialCorrect    = []
    # trialLabel      = []
    # trialFeedback   = []
    #
    # # Trial evaluation
    # # =========================================================================
    #
    # ix = log.index.tolist()[0]
    #
    # # Match pattern, using stimulus and response data
    # source = eval_data['eval_data'].fillna(float('inf'))
    # patternDict= log[source.columns].fillna(float('inf')).to_dict()
    # pattern = {key: [value[ix]] for key,value in patternDict.iteritems()}
    # iRow = source.isin(pattern).all(1)
    #
    # if sum(iRow) == 0:
    #
    #     trialCorrect = False
    #
    #     # If no match, try to match using stimulus data only to determine trialType
    #     sourceStim = eval_data['eval_data'][['s1Ix','s2Ix']].fillna(float('inf'))
    #     patternDictStim = log[sourceStim.columns].fillna(float('inf')).to_dict()
    #     patternStim = {key: [value[ix]] for key,value in patternDictStim.iteritems()}
    #     iRow = sourceStim.isin(patternStim).all(1)
    #
    #     uniqueTrialTypes = eval_data['trialType'].loc[iRow].unique()
    #     if len(uniqueTrialTypes) == 0:
    #         trialType = 'NOC' # Not otherwise classified
    #     elif len(uniqueTrialTypes) == 1:
    #         trialType = uniqueTrialTypes.tolist()[0]
    #     else:
    #         trialType = "Non-unique: %s" % ', '.join(uniqueVals.tolist())
    #
    #     responseType = 'NOC' # Not otherwise classified
    #     trialFeedback = 'incorrect'
    #
    # elif sum(iRow) == 1:
    #
    #     trialCorrect    = eval_data['correct'].loc[iRow].as_matrix().tolist()[0]
    #     trialType       = eval_data['trialType'].loc[iRow].as_matrix().tolist()[0]
    #     responseType    = eval_data['responseType'].loc[iRow].as_matrix().tolist()[0]
    #     trialFeedback   = eval_data['feedback'].loc[iRow].as_matrix().tolist()[0]
    #
    # elif sum(iRow) > 1:
    #
    #     # Determine if there is one unique trialCorrect value among selected rows
    #     uniqueTrialCorrect = eval_data['correct'].loc[iRow].unique()
    #     if len(uniqueTrialCorrect) == 0:
    #         trialCorrect = 'Unknown'
    #     elif len(uniqueTrialCorrect) == 1:
    #         trialCorrect = uniqueTrialCorrect.tolist()[0]
    #     else:
    #         trialCorrect = "Non-unique: %s" % ', '.join(uniqueTrialCorrect.tolist())
    #
    #     # Determine if there is one unique trialType among selected rows
    #     uniqueTrialTypes = eval_data['trialType'].loc[iRow].unique()
    #     if len(uniqueTrialTypes) == 0:
    #         trialType = 'NOC' # Not otherwise classified
    #     elif len(uniqueTrialTypes) == 1:
    #         trialType = uniqueTrialTypes.tolist()[0]
    #     else:
    #         trialType = "Non-unique: %s" % ', '.join(uniqueTrialTypes.tolist())
    #
    #     # Determine if there is one unique response type among selected rows
    #     uniqueResponseTypes = eval_data['responseType'].loc[iRow].unique()
    #     if len(uniqueResponseTypes) == 0:
    #         responseType = 'NOC' # Not otherwise classified
    #     elif len(uniqueResponseTypes) == 1:
    #         responseType = uniqueResponseTypes.tolist()[0]
    #     else:
    #         responseType = "Non-unique: %s" % ', '.join(uniqueResponseTypes.tolist())
    #
    #     # Determine if there is one unique trial feedback among selected rows
    #     uniqueTrialFeedback = eval_data['feedback'].loc[iRow].unique()
    #     if len(uniqueTrialFeedback) == 0:
    #         trialFeedback = 'incorrect' # Not otherwise classified
    #     elif len(uniqueTrialFeedback) == 1:
    #         trialFeedback = uniqueTrialFeedback.tolist()[0]
    #     else:
    #         trialFeedback = 'incorrect'
    #
    #     print('Warning: combination of stimuli and response(s) matches multiple criteria.')
    #     print('--------------------------------------------------------------')
    #     print('Stimuli and response(s) on this trial:')
    #     pprint(pattern)
    #     print('trialCorrect: \t %s' % (trialCorrect))
    #     print('trialType: \t %s' % (trialType))
    #     print('responseType: \t %s' % (responseType))
    #     print('trialFeedback: \t %s' % (trialFeedback))
    #
    # log.iloc[0]['trialCorrect']     = trialCorrect
    # log.iloc[0]['trialType']        = trialType
    # log.iloc[0]['responseType']     = responseType
    # log.iloc[0]['trialFeedback']    = trialFeedback

    # SOA adjustments
    # =========================================================================

    # newSoa = []
    #
    # if soa['type'] == 'staircase':
    #
    #     soaIx = []
    #     soa = []
    #
    #
    #     if trialLabel == 'correct':
    #         newSoa = soa + soaStep
    #     else:
    #         newSoa = soa - soaStep

    # Feedback
    # =========================================================================

    # stimuli['feedback'][0].setText(trialFeedback)
    # stimuli['feedback'][0].setAutoDraw(True)
    # stimuli['iti'][0].setAutoDraw(True)
    # window.flip()
    # pp.core.wait(feedback_dur)
    # stimuli['feedback'][0].setAutoDraw(False)
    stimuli['iti'][0].setAutoDraw(True)
    stimuli['iti'][0].setText('o')
    window.flip()
    stimuli['iti'][0].setAutoDraw(False)


    return log
def init_config(runtime, mod_dir):
    """
    Parse and process experiment configuration file.

    For details about its contents, see the experiment configurations file.

    Parameters
    ----------
    runtime     : __main__.ExperimentRuntime
                ioHubExperimentRuntime object, which brings together several
                aspects of the ioHub Event Monitoring Framework.

    mod_dir     : str or unicode
                directory where itch_task module resides.

    Returns
    -------
    config      : dict
                specifies itch_task experiment properties

    """

    config = runtime.getConfiguration()

    user_def_params = runtime.getUserDefinedParameters()

    hub = runtime.hub

    ###########################################################################
    # STUDY
    ###########################################################################

    config['study'] = dict(study_id=config['code'],
                           task_version_id=config['version'],
                           ethics_protocol_id=config['ethics_protocol_id'])

    ###########################################################################
    # SUBJECT
    ###########################################################################
    # subject_ix
    # group_ix
    # sex
    # age

    config['subject'] = {'subject_ix': int(user_def_params['subject_ix']),
                         'group_ix': int(user_def_params['group_ix'])}

    ###########################################################################
    # SESSION
    ###########################################################################

    t_sess_start_unix = time.time() - pp.core.getTime()  # UNIX time stamp
    t_sess_start_gmt = time.gmtime(t_sess_start_unix)  # GMT time stamp

    # Seed the random number generator
    random.seed(t_sess_start_unix)

    session = {'session_ix': int(user_def_params['session_ix']),
               'session_id': config['session_defaults']['name'],
               'experimenter_id': config['session_defaults'][
                   'experimenter_id'],
               'date': time.strftime('%Y-%m-%d', t_sess_start_gmt),
               'weekday': time.strftime('%a', t_sess_start_gmt),
               'time': time.strftime('%H%M-GMT', t_sess_start_gmt),
               'rng_seed': t_sess_start_unix}

    config['session'] = session

    ###########################################################################
    # LOG - PART 1: CHECK IF THE EXPERIMENTER ENTERED THE CORRECT DATA
    ###########################################################################

    group_ix = config['subject']['group_ix']
    subject_ix = config['subject']['subject_ix']
    session_ix = config['session']['session_ix']
    study_id = config['study']['study_id']
    task_version_id = config['study']['task_version_id']
    expt_dir = os.path.abspath(os.path.join(mod_dir, os.pardir))

    # Make a log directory, if it does not exist
    if config['log']['dir']:
        log_dir = os.path.normcase(
            os.path.join(config['log']['dir'], study_id))
    else:
        log_dir = os.path.normcase(os.path.join(expt_dir, 'log/', study_id))

    if not os.path.isdir(log_dir):
        os.makedirs(log_dir)

    str_format_performance = \
        '%s_Study_%s_TaskVersion_%s_Group_%.2d_Subject_%.3d'

    trial_log_file_name = str_format_performance % ('trial_log',
                                                    study_id,
                                                    task_version_id,
                                                    group_ix,
                                                    subject_ix)
    trial_log_file = \
        os.path.normcase(
            os.path.join(log_dir,
                         trial_log_file_name + '.csv')
        )

    # If the file exists, check if it contains data corresponding to session_ix
    if os.path.isfile(trial_log_file):
        group_ix = config['subject']['group_ix']
        subject_ix = config['subject']['subject_ix']

        trial_log = pd.read_csv(trial_log_file,
                                error_bad_lines=False)

        if session_ix in trial_log.session_ix.values:
            warn_dlg = pp.gui.Dlg(title="WARNING",
                                  labelButtonOK=u' Continue ',
                                  labelButtonCancel=u' Cancel ')
            warn_dlg.addText('You specified the following settings:')
            warn_dlg.addFixedField('Group index:', group_ix)
            warn_dlg.addFixedField('Subject index:', subject_ix)
            warn_dlg.addFixedField('Session index:', session_ix)
            warn_dlg.addText('')
            # TODO: Check how this text can be wrapped
            warn_dlg.addText(
                'You might have entered the wrong data. A log file with these data already exists:')
            warn_dlg.addText(trial_log_file)
            warn_dlg.addText('')
            # TODO: Check how this text can be wrapped
            warn_dlg.addText(
                'Press Continue if you want to use the above settings and overwrite/append this file.')
            warn_dlg.addText('Press Cancel if you want to change settings.')

            warn_dlg.show()

            if not warn_dlg.OK:
                return -1

    ###########################################################################
    # APPARATUS
    ###########################################################################

    config['apparatus'] = {'hub': hub,
                           'display': dict(),
                           'kb': dict(),
                           'rd': dict()}

    config['apparatus']['display']['client'] = hub.getDevice('display')
    config['apparatus']['kb']['client'] = hub.getDevice('keyboard')

    if hub.getDevice('response_device') is None:
        config['apparatus']['rd']['client'] = hub.getDevice('keyboard')
    else:
        config['apparatus']['rd']['client'] = hub.getDevice('response_device')

    # Keyboard settings
    # -------------------------------------------------------------------------
    esc_keys = config['responses']['abort_keys']
    trigger_keys = config['responses']['trigger_keys']
    toggle_keys = config['responses']['toggle_keys']
    config['apparatus']['kb']['settings'] = {'esc_keys': esc_keys,
                                             'trigger_keys': trigger_keys,
                                             'toggle_keys': toggle_keys}

    kb = config['apparatus']['kb']['client']

    # Response device settings
    # -------------------------------------------------------------------------
    rd = config['apparatus']['rd']['client']
    rd_class = rd.getIOHubDeviceClass()
    rsp_keys_per_class = config['responses']['response_keys_per_class']
    key_key = {'Keyboard': 'key',
               'Serial': 'data'}
    time_key = {'Keyboard': 'time',
                'Serial': 'time'}

    config['apparatus']['rd']['settings'] = {'class': rd_class,
                                             'rsp_keys': rsp_keys_per_class[
                                                 rd_class],
                                             'key_key': key_key[rd_class],
                                             'time_key': time_key[rd_class]}

    # Enable event reporting and clear all recorded events
    rd.enableEventReporting()
    rd.clearEvents()

    ###########################################################################
    # WINDOW
    ###########################################################################
    display = config['apparatus']['display']['client']

    if __debug__:
        window = pp.visual.Window(display.getPixelResolution(),
                                  monitor=display.getPsychopyMonitorName(),
                                  units='deg',
                                  fullscr=False,
                                  allowGUI=False,
                                  screen=display.getIndex())
    else:
        window = pp.visual.Window(display.getPixelResolution(),
                                  monitor=display.getPsychopyMonitorName(),
                                  units='deg',
                                  fullscr=True,
                                  allowGUI=False,
                                  screen=display.getIndex())

    frame_rate = None
    cntr = 0
    while frame_rate is None:
        cntr = cntr + 1
        frame_rate = window.getActualFrameRate()

        if cntr >= 3:
            print 'Cannot get frame rate after %d attempts.' % cntr
            print 'Terminating PsychoPy. Try and restart the program.'
            pp.core.quit()

    frame_time = 1 / frame_rate
    config['window'] = {'window': window,
                        'frame_rate': frame_rate,
                        'frame_time': frame_time}

    ###########################################################################
    # STIMULI
    ###########################################################################

    # Make stimuli
    config['stimuli'] = {
        stim: define_stimulus(window, config['stim_config'][stim], stim)
        for stim in config['stim_config']}

    ###########################################################################
    # RESPONSES
    ###########################################################################

    ###########################################################################
    # EVALUATION
    ###########################################################################
    # TODO: Re-implement evaluation

    # trial_eval_data = pd.read_csv(
    #     config['evaluation']['trial']['eval_data_file'][rd_class],
    #     error_bad_lines=False)
    # trial_eval_data = check_df_from_csv_file(trial_eval_data)
    #
    # trial_categories = trial_eval_data.fillna(value=np.nan)
    #
    # eval_columns = []
    # for col in trial_eval_data.columns:
    #     if not col.startswith('trial'):
    #         eval_columns.append(col)
    #
    # config['evaluation']['trial']['eval_data'] = trial_eval_data[
    #     eval_columns].copy()
    # config['evaluation']['trial']['correct'] = trial_eval_data[
    #     'trial_correct'].copy()
    # config['evaluation']['trial']['trial_type'] = trial_eval_data[
    #     'trial_type'].copy()
    # config['evaluation']['trial']['response_type'] = trial_eval_data[
    #     'trialresponse_type'].copy()
    # config['evaluation']['trial']['feedback'] = trial_eval_data[
    #     'trial_feedback'].copy()

    ###########################################################################
    # INSTRUCTION
    ###########################################################################

    # config['stimuli']['instruction'] = \
    #     {instr_type: define_stimulus(window,
    #                                  config['instruction'][instr_type])
    #      for instr_type in config['instruction']
    #      if config['instruction'][instr_type]['content'] is not None}
    #
    # instr_list_p = pd.read_csv(
    #     config['instruction']['practice']['instruction_list_file'],
    #     error_bad_lines=False)
    # instr_list_p = check_df_from_csv_file(df=instr_list_p)
    # config['instruction']['practice']['list'] = instr_list_p
    #
    # instr_list_e = pd.read_csv(
    #     config['instruction']['experiment']['instruction_list_file'],
    #     error_bad_lines=False)
    # instr_list_e = check_df_from_csv_file(df=instr_list_e)
    # config['instruction']['experiment']['list'] = instr_list_e

    ###########################################################################
    # PRACTICE
    ###########################################################################

    if user_def_params['practice']:
        config['practice']['enable'] = True
    else:
        config['practice']['enable'] = False

    # If a trial_listFile exists, use this
    tr_list_p = pd.read_csv(config['practice']['trial_list_file'],
                          error_bad_lines=False)
    # tr_list_p = check_df_from_csv_file(df=tr_list_p)
    config['practice']['trial_list'] = tr_list_p

    ###########################################################################
    # EXPERIMENT
    ###########################################################################

    if user_def_params['experiment']:
        config['experiment']['enable'] = True
    else:
        config['experiment']['enable'] = False

    tr_list_e = pd.read_csv(config['experiment']['trial_list_file'],
                            error_bad_lines=False)
    # tr_list_e = check_df_from_csv_file(df=tr_list_e)
    config['experiment']['trial_list'] = tr_list_e

    ###########################################################################
    # INDIFFERENCE POINT PROCEDURE
    ###########################################################################

    ###########################################################################
    # PERFORMANCE REQUIREMENTS
    ###########################################################################

    ###########################################################################
    # CLOCK
    ###########################################################################

    config['clock'] = pp.core.Clock()

    ###########################################################################
    # LOG - PART 2: INITIATE LOG FILES
    ###########################################################################

    # Run time info
    # -------------------------------------------------------------------------
    str_format_runtime = \
        '%s_Study_%s_Group_%.2d_Subject_%.3d_Session_%.2d_%s_%s'

    if config['log']['runtime']['enable']:
        run_time_info = ppinfo.RunTimeInfo(win=window,
                                           refreshTest=True,
                                           verbose=True,
                                           userProcsDetailed=True)

        run_time_info_file_name = str_format_runtime % ('run_time_info',
                                                    config['study'][
                                                        'study_id'],
                                                    config['subject'][
                                                        'group_ix'],
                                                    config['subject'][
                                                        'subject_ix'],
                                                    config['session'][
                                                        'session_ix'],
                                                    config['session']['date'],
                                                    config['session']['time'])

        run_time_info_file = os.path.normcase(
            os.path.join(log_dir, run_time_info_file_name + '.csv'))

        with open(run_time_info_file, 'a+') as fileObj:
            fileObj.write(str(run_time_info))

    # Task performance
    # -------------------------------------------------------------------------

    # Init log data frame
    trial_cols, block_cols, sess_cols, sess_data = init_log(config)

    config['log']['performance']['sess_columns'] = sess_cols
    config['log']['performance']['sess_data'] = sess_data

    if config['log']['performance']['trial']['enable']:

        config['log']['performance']['trial']['columns'] = trial_cols

        trial_log_file_name = str_format_performance % ('trial_log',
                                                        config['study'][
                                                            'study_id'],
                                                        config['study'][
                                                            'task_version_id'],
                                                        config['subject'][
                                                            'group_ix'],
                                                        config['subject'][
                                                            'subject_ix'])

        trial_log_file = os.path.normcase(
            os.path.join(log_dir, trial_log_file_name + '.csv'))

        config['log']['performance']['trial']['file'] = trial_log_file

        if not os.path.isfile(trial_log_file):
            with open(trial_log_file, 'a+') as fileObj:
                pd.DataFrame(index=[], columns=trial_cols).to_csv(fileObj,
                                                                  index=False,
                                                                  header=True)

    if config['log']['performance']['block']['enable']:

        config['log']['performance']['block']['columns'] = block_cols

        block_log_file_name = str_format_performance % ('block_log',
                                                        config['study'][
                                                            'study_id'],
                                                        config['study'][
                                                            'task_version_id'],
                                                        config['subject'][
                                                            'group_ix'],
                                                        config['subject'][
                                                            'subject_ix'])

        block_log_file = os.path.normcase(
            os.path.join(log_dir, block_log_file_name + '.csv'))

        config['log']['performance']['block']['file'] = block_log_file

        if not os.path.isfile(block_log_file):
            with open(block_log_file, 'a+') as fileObj:
                pd.DataFrame(index=[], columns=block_cols).to_csv(fileObj,
                                                                  index=False,
                                                                  header=True)

    return config
def init_log(config):
    """
    Initializes and saves a pandas DataFrame for logging stop task data and
    returns

    Parameters
    ----------
    config      : dict
                specifies StPy experiment properties

    Returns
    -------
    trial_cols  : list
                columns in trial log

    block_cols  : list
                columns in block log

    sess_columns : list
                columns in trial log and block log containing session-specific
                data

    sess_data   : list
                session-specific data


    """

    # Process inputs
    # -------------------------------------------------------------------------

    # Assertions

    ###########################################################################
    # SESSION-SPECIFIC COLUMNS AND DATA
    ###########################################################################

    id_columns_sess = ['study_id',
                       'task_version_id',
                       'session_id',
                       'experimenter_id']

    id_data_sess = [config['study']['study_id'],
                    config['study']['task_version_id'],
                    config['session']['session_id'],
                    config['session']['experimenter_id']]

    tech_columns_sess = ['response_device',
                         'refresh_rate',
                         'rng_seed']

    tech_data_sess = [config['apparatus']['rd']['settings']['class'],
                      config['window']['frame_rate'],
                      config['session']['rng_seed']]

    ix_columns_sess = ['subject_ix',
                       'group_ix',
                       'session_ix']

    ix_data_sess = [config['subject']['subject_ix'],
                    config['subject']['group_ix'],
                    config['session']['session_ix']]

    time_columns_sess = ['sess_date',
                         'sess_time']
    time_data_sess = [config['session']['date'],
                      config['session']['time']]

    sess_columns = []
    sess_columns += id_columns_sess
    sess_columns += tech_columns_sess
    sess_columns += ix_columns_sess
    sess_columns += time_columns_sess

    sess_data = []
    sess_data += id_data_sess
    sess_data += tech_data_sess
    sess_data += ix_data_sess
    sess_data += time_data_sess

    # # Copy columns to prevent that sess_columns is updated
    # trialColumns    = list(sess_columns)
    # blockColumns    = list(sess_columns)

    ###########################################################################
    # TRIAL-SPECIFIC COLUMNS
    ###########################################################################

    def init_trial_log(config, columns):

        # Identifiers
        # =====================================================================
        id_columns = ['block_id']

        # Indices
        # =====================================================================
        ix_columns = ['block_ix',
                      'trial_ix']

        # Framing
        # =====================================================================
        framing_columns = ['framing']

        # Position
        # =====================================================================
        pos_columns = ['ll_side']

        # Time
        # =====================================================================
        time_columns = ['t_session',
                        't_block',
                        'trial_ons',
                        'trial_dur',
                        'choice_instr_ons',
                        'choice_instr_ons_dt',
                        'choice_instr_dur',
                        'choice_instr_dur_dt',
                        'ss_opt_ons',
                        'ss_opt_ons_dt',
                        'ss_opt_dur',
                        'ss_opt_dur_dt',
                        'll_opt_ons',
                        'll_opt_ons_dt',
                        'll_opt_dur',
                        'll_opt_dur_dt']

        # Trigger
        # =====================================================================
        trigger_columns = ['waited_for_trigger']

        # Feedback
        # =====================================================================
        # TODO: Think about feedback logging
        # feedback_columns = ['trial_correct',
        #                     'trial_type',
        #                     'response_type',
        #                     'trial_feedback']

        # Choice options
        # =====================================================================
        choice_opt_columns = ['trial_type', 'm_s', 'm_l', 'm_unit', 't_s',
                              't_l', 't_unit']

        # Response event times
        # =====================================================================
        rsp_keys = [item
                    for sublist
                    in config['apparatus']['rd']['settings']['rsp_keys']
                    for item in sublist]

        resp_ev_columns = ['key_count_' + key for key in rsp_keys]
        resp_ev_columns += ['choice',
                            'rt'
                            ]

        # Put it together
        # =====================================================================
        columns += id_columns
        columns += ix_columns
        columns += framing_columns
        columns += pos_columns
        columns += time_columns
        columns += trigger_columns
        columns += choice_opt_columns
        columns += resp_ev_columns
        # columns += feedback_columns

        return columns

    trial_cols = init_trial_log(config=config, columns=list(sess_columns))

    ###########################################################################
    # BLOCK-SPECIFIC COLUMNS
    ###########################################################################

    def init_block_log(config, columns):

        # Identifiers
        # =====================================================================
        id_columns = ['block_id']

        # Indices
        # =====================================================================
        ix_columns = ['block_ix',
                      'iterIx']

        # n_s1 = len(config['stimuli']['s1'])
        # n_s2 = len(config['stimuli']['s2'])

        # Accuracy - s1
        # =====================================================================
        # if config['feedback']['block']['features']['s1_accuracy']['enable']:
        #
        #     s1_acc_comprehension = [('s1_acc_%.2d' % s1,
        #                              's1_acc_crit_met_%.2d' % s1) for s1 in
        #                             range(n_s1)]
        #     s1_acc_cols = list(it.chain(*s1_acc_comprehension))
        #
        # else:
        #     s1_acc_cols = []

        # Accuracy - s2
        # =====================================================================
        # if config['feedback']['block']['features']['s2_accuracy']['enable']:
        #     s2_acc_comprehension = [('s2_acc_%.2d' % s2,
        #                              's2_acc_crit_met_%.2d' % s2) for s2 in
        #                             range(n_s2)]
        #     s2_acc_cols = list(it.chain(*s2_acc_comprehension))
        # else:
        #     s2_acc_cols = []

        # Mean RT - s1
        # =====================================================================
        # if config['feedback']['block']['features']['s1MeanRt']['enable']:
        #     s1_rt_comprehension = [('s1_mean_rt_%.2d' % s1,
        #                             's1_mean_rt_crit_met_%.2d' % s1) for s1 in
        #                            range(n_s1)]
        #     s1_rt_cols = list(it.chain(*s1_rt_comprehension))
        # else:
        #     s1_rt_cols = []

        # Mean RT diff - s1
        # =====================================================================
        # if config['feedback']['block']['features']['s1_meanrt_diff']['enable']:
        #     s1_rt_diff_comprehension = \
        #         [('s1_mean_rt_diff_%.2d' % s1,
        #           's1_mean_rt_diff_crit_met_%.2d' % s1)
        #          for s1 in range(n_s1)]
        #     s1_rt_diff_cols = list(it.chain(*s1_rt_diff_comprehension))
        # else:
        #     s1_rt_diff_cols = []

        # Put it together
        # =====================================================================
        columns += id_columns
        columns += ix_columns
        # columns += s1_acc_cols
        # columns += s2_acc_cols
        # columns += s1_rt_cols
        # columns += s1_rt_diff_cols

        return columns

    block_cols = init_block_log(config=config, columns=list(sess_columns))

    ###########################################################################
    # OUTPUT
    ###########################################################################

    return trial_cols, block_cols, sess_columns, sess_data
def init_stimulus(window, stim_type):
    """

    Initialize PsychoPy stimulus with default settings

    Parameters
    ----------
    window      : psychopy.visual.window.Window
                PsychoPy window object, in which stimuli are presented

    stim_type    : str or unicode
                stimulus type
    Returns
    -------
    stim_object  : psychopy.visual.text.TextStim or psychopy.visual.text.ImageStim
                PsychoPy stimulus object
    """

    stim_dict = {'textstim': pp.visual.TextStim(window),
                 'imagestim': pp.visual.ImageStim(window)}

    stim_object = stim_dict[stim_type.lower()]

    return stim_object
def present_instruction(config, block_type, *args):
    """
    Present instruction stimuli

    Parameters
    ----------
    config              : dict
                        specifies StPy experiment properties

    block_type          : str or unicode
                        block type

    block_ix            : int (optional)
                        block index

    instruction_stim_ix : int (optional)
                        instruction stimulus index

    """

    # Process inputs
    # -------------------------------------------------------------------------

    # Process variable arguments
    if len(args) > 0:
        if block_type == 'practice' or block_type == 'experiment':
            block_ix = args[0]
        elif block_type == 'block_repeat':
            instruction_stim_ix = args[0]
        else:
            block_ix = None
            instruction_stim_ix = None

    window = config['window']['window']
    hub = config['apparatus']['hub']
    rd = config['apparatus']['rd']
    kb = config['apparatus']['kb']
    toggle_keys = kb['settings']['toggle_keys']
    instruction_stim = config['stimuli']['instruction'][type]

    stim_ix = 0

    # Determine which instruction screens to show, depending on the
    # experiment phase
    if type == 'practice' or type == 'experiment':
        session_ix = config['session']['session_ix']
        instruction_list = config['instruction'][type]['list']

        pattern = {'session_ix': [session_ix],
                   'block_ix': [block_ix]}
        ix = instruction_list.index.tolist()
        i_row = instruction_list[pattern.keys()].isin(pattern).all(1)
        stim_list = instruction_list.loc[i_row, 'instruction_ix'].astype(
            int).tolist()
    elif type == 'block_repeat':
        stim_list = [instruction_stim_ix]
    elif type == 'start' or type == 'end':
        stim_list = range(len(instruction_stim))

    # Show instruction screens and monitor keyboard and response device inputs
    # -------------------------------------------------------------------------

    while stim_ix < len(stim_list):

        # Clear buffer
        hub.clearEvents('all')
        rd['client'].clearEvents()
        kb['client'].clearEvents()
        no_key_pressed = True

        # Show the instruction
        instruction_stim[stim_list[stim_ix]].draw()
        window.flip()

        while no_key_pressed:

            # Collect responses
            rd_key_count, toggle_keys_pressed = \
                collect_response(rd, kb, other_keys=toggle_keys)

            # If user pressed key move to next stimulus
            if sum(rd_key_count.values()) > 0:
                window.flip(clearBuffer=True)
                stim_ix += 1
                break

            # If toggle keys are used move to next or previous stimulus
            if toggle_keys_pressed:
                window.flip(clearBuffer=True)
                if toggle_keys_pressed == toggle_keys[0]:
                    stim_ix -= 1
                    if stim_ix < 0:
                        stim_ix = 0
                    break
                elif toggle_keys_pressed == toggle_keys[1]:
                    stim_ix += 1
                    break
def present_stimuli(window, stim_list, u, f_on_off, log):
    """
    Presents stimuli comprising a trial (choice_instr, ss_opt, and ll_opt)

    Parameters
    ----------
    window      : psychopy.visual.window.Window
                PsychoPy window object, in which stimuli are presented

    stim_list   : list
                specifies a list of PsychoPy stimuli to present on this trial

    u           : numpy.ndarray
                stimulus-by-frame array specifying for each combination of
                stimulus (row) and frame (column) whether or not the stimulus
                should be displayed

    f_on_off    : numpy.ndarray
                stimulus-by-2 array specifying the frame indices of stimulus
                onset (first column) and stimulus offset (second column)

    log         : pandas.core.frame.DataFrame
                trial log

    Returns
    -------
    log         : pandas.core.frame.DataFrame
                trial log; present_stimuli fills in information about trial and
                stimulus onsets and durations, as well as differences between
                actual and planned onsets and durations

    """

    # Define dynamic variables
    # -------------------------------------------------------------------------
    n_stim = len(stim_list)
    n_frame = np.size(u, 1)
    t_flip = [None] * (n_frame + 1)

    # Execute trial
    # =========================================================================

    # Draw stimuli, flip screen
    for frame_ix in range(n_frame):

        for stim_ix in range(n_stim - 1):

            if u[stim_ix, frame_ix]:
                stim_list[stim_ix].setAutoDraw(True)
            else:
                stim_list[stim_ix].setAutoDraw(False)

        t_flip[frame_ix] = window.flip()

    # Hide all trial stimuli and present ITI stimulus
    [stim_list[stim_ix].setAutoDraw(False) for stim_ix in range(n_stim - 1)]
    stim_list[-1].setAutoDraw(True)
    t_flip[-1] = window.flip()

    # Timing
    trial_ons = t_flip[0]
    trial_off = t_flip[-1]
    trial_dur = trial_off - trial_ons

    # Actual stimulus onset and duration times
    stim_displayed = [stim.name for stim in stim_list]

    # Log
    log.iloc[0]['trial_ons'] = trial_ons
    log.iloc[0]['trial_dur'] = trial_dur

    for ix in range(len(stim_displayed) - 1):
        f_on, f_off = f_on_off[:, ix]
        ons = t_flip[f_on] - trial_ons
        dur = t_flip[f_off] - t_flip[f_on]
        ons_intended = log.iloc[0][stim_displayed[ix] + '_ons']
        dur_intended = log.iloc[0][stim_displayed[ix] + '_dur']
        ons_dt = ons - ons_intended
        dur_dt = dur - dur_intended

        log.iloc[0][stim_displayed[ix] + '_ons'] = ons
        log.iloc[0][stim_displayed[ix] + '_dur'] = dur
        log.iloc[0][stim_displayed[ix] + '_ons_dt'] = ons_dt
        log.iloc[0][stim_displayed[ix] + '_dur_dt'] = dur_dt

    return log

def run_block(config,block_id,trial_list,block_log,trial_ons_next_block):
    """
    Runs all trials comprising a block

    Parameters
    ----------
    config              : dict
                        specifies StPy experiment properties

    block_id             : str or unicode
                        identifier of the block

    trial_list           : pandas.core.frame.DataFrame
                        list of trials to present

    block_log            : pandas.core.frame.DataFrame
                        log of trials comprising the present block

    trial_ons_next_block  : float
                        onset of the next block (in seconds, relative to last
                        clock reset)


    Returns
    -------
    block_log            : pandas.core.frame.DataFrame
                        block-level performance log

    all_crit_met          : bool
                        whether or not all performance criteria have been met

    """

    # 1. Process inputs
    # =========================================================================

    # All relevant variables and objects
    window  = config['window']['window']
    stimuli = config['stimuli']
    hub = config['apparatus']['hub']
    rd = config['apparatus']['rd']
    kb = config['apparatus']['kb']
    trial_stats = config['statistics']['trial']
    ip_procedure = config['ip_procedure']

    trial_eval_data = config['evaluation']['trial']
    feedback_dur = config['feedback']['trial']['duration']
    sess_columns = config['log']['performance']['sess_columns']
    sess_data = config['log']['performance']['sess_data']

    block_ix = trial_list.iloc[0]['block_ix']


    # If this is the indifference point procedure
    # if block_id.startswith('i'):
    #     m_l = float(config['ip_procedure']['m_l'])
    #     m_s = m_l * 0.5
    #     m_units = config['ip_procedure']['m_units']
    #     t_l_list = config['ip_procedure']['t_l']
    #     t_s = int(config['ip_procedure']['t_s'])
    #     t_units = config['ip_procedure']['t_units']

    # Define dynamic variables
    # -------------------------------------------------------------------------

    # Mock times, so that trial starts immediately
    trial_timing = {'ons': -float('inf'),
                    'dur': 0,
                    'iti_dur': 0,
                    'refresh_time': 1/config['window']['frame_rate']}
    block_ons = -float('inf')

    # =========================================================================
    trial_list_ixs = trial_list.index.tolist()
    trial_cols = config['log']['performance']['trial']['columns']
    trial_log = pd.DataFrame(index = trial_list_ixs,
                             columns=trial_cols)

    # If this is the indifference point procedure
    if block_id.startswith('i'):

        i_step = iter(range(ip_procedure['n_staircase_trial']))
        m_l = ip_procedure['m_l']

        step_size = m_l * 2 ** -(next(i_step) + 1)
        adjustment_factor = -1
        m_s = m_l + adjustment_factor * step_size

    # Present trials
    # =========================================================================

    for trial_list_ix in trial_list_ixs:

        # Prepare trial
        # ---------------------------------------------------------------------
        # TODO: Change to pd.Series? e.g. this_trial_log = pd.Series(index = trial_cols)
        this_trial_log = pd.DataFrame(index = [trial_list_ix],
                                      columns = trial_cols)

        # If this is the indifference point procedure
        if block_id.startswith('i'):
            if trial_list.loc[trial_list_ix, 'trial_type'] in ['catch_ss',
                                                               'catch_ll',
                                                               'instr_check']:

                # Set monetary amount of SS option
                m_s = {'catch_ss': ip_procedure['m_l'],
                       'catch_ll': 0,
                       'instr_check': m_s}

            elif trial_list.loc[trial_list_ix, 'trial_type'] == 'standard':
                # TODO: Set m_s to correct amount
                m_s = m_s

        this_trial_log.iloc[0]['trial_type'] = \
            trial_list.loc[trial_list_ix, 'trial_type']

        # Check if program should wait for trigger
        if 'wait_for_trigger' in trial_list.columns:
            wait_for_trigger = trial_list.ix[trial_list_ix]['wait_for_trigger']
        else:
            wait_for_trigger = False

        trial_ons = trial_list.loc[trial_list_ix,'trial_ons']

        # Randomly determine side of SS and LL option
        x_ss, y_ss = config['stimuli']['ss_opt'][0].pos
        x_ll, y_ll = config['stimuli']['ll_opt'][0].pos

        if random.random() < 0.5:
            stimuli['ss_opt'][0].pos = (x_ss, y_ss)
            stimuli['ll_opt'][0].pos = (x_ll, y_ll)
            this_trial_log.iloc[0]['ll_side'] = 'right'
        else:
            # Swap sides
            stimuli['ss_opt'][0].pos = (x_ll, y_ll)
            stimuli['ll_opt'][0].pos = (x_ss, y_ss)
            this_trial_log.iloc[0]['ll_side'] = 'left'

        for stimulus in ['choice_instr', 'ss_opt', 'll_opt']:
            stimuli[stimulus][0].alignHoriz = 'center'
            stimuli[stimulus][0].alignVert = 'center'
            stimuli[stimulus][0].wrapWidth = 1000

            stimuli[stimulus][0].setText(
                make_stim_str(stim=stimulus,
                              framing=trial_list.loc[trial_list_ix,'framing'],
                              m_unit=trial_list.loc[trial_list_ix,'m_units'],
                              m_s=m_s,
                              m_l=trial_list.loc[trial_list_ix,'m_l'],
                              t_unit=trial_list.loc[trial_list_ix,'t_units'],
                              t_s=trial_list.loc[trial_list_ix,'t_s'],
                              t_l=trial_list.loc[trial_list_ix,'t_l'])
            )

        # Fill in session data
        # ---------------------------------------------------------------------
        this_trial_log.loc[trial_list_ix,sess_columns] = sess_data

        this_trial_log, stim_list, u, f_on_off, trial_dur = \
            stim_to_frame_mat(config,
                              trial_list.ix[trial_list_ix],
                              this_trial_log)

        t_trial_ready = config['clock'].getTime()
        print 'Trial %d ready to start: t = %f s, dt = %f ms' % (trial_list_ix, t_trial_ready, 1000*(t_trial_ready - trial_ons))

        # Run trial
        # ---------------------------------------------------------------------
        this_trial_log = run_trial(config,
                                   wait_for_trigger,
                                   trial_ons,
                                   hub,
                                   this_trial_log,
                                   trial_timing,
                                   window,
                                   stim_list,
                                   u,
                                   f_on_off,
                                   rd,
                                   kb,
                                   trial_stats,
                                   trial_eval_data,
                                   feedback_dur,
                                   stimuli)

        # Log trial data not logged inside run_trial
        # ---------------------------------------------------------------------
        this_trial_log.loc[trial_list_ix,'block_id'] = block_id
        this_trial_log.loc[trial_list_ix,'block_ix'] = block_ix
        this_trial_log.loc[trial_list_ix,'trial_ix'] = trial_list.ix[
            trial_list_ix]['trial_ix']

        # Session timing
        sm, ss = divmod(this_trial_log['trial_ons'].item(), 60)
        sh, sm = divmod(sm, 60)
        this_trial_log.loc[trial_list_ix,'t_session']    = '%d:%02d:%02d' % (sh, sm, ss)

        # Block timing
        if trial_list_ix == trial_list_ixs[0]:
            block_ons = this_trial_log['trial_ons'].item()
        bs, bms = divmod(this_trial_log['trial_ons'].item() - block_ons,1)
        bm, bs = divmod(bs, 60)
        bh, bm = divmod(bm, 60)
        this_trial_log.loc[trial_list_ix,'tBlock']      = '%d:%02d:%02d.%03d' % (bh, bm, bs,bms*1000)

        # Put trial data into data frame and file
        # ---------------------------------------------------------------------
        trial_log.iloc[trial_list_ix] = this_trial_log.iloc[trial_list_ix]

        trial_ons = time.time();

        with open(config['log']['performance']['trial']['file'],'a+') as fileObj:
            this_trial_log.to_csv(fileObj, index=False, header=False, na_rep=np.nan)

        print 'Time needed to write trial_log: %.3f ms' % (1000.*(time.time() - trial_ons))

        # Update m_s for upcoming trial
        # ---------------------------------------------------------------------
        if block_id.startswith('i'):

            if trial_list.loc[trial_list_ix, 'trial_type'] == 'standard':

                step_size = m_l * 2 ** -(next(i_step) + 1)

                adjustment_factor = \
                    {'ss': -1,
                     'll': 1}[this_trial_log.loc[trial_list_ix, 'choice']]

                m_s = m_s + adjustment_factor * step_size


    # Compute block stats
    # =========================================================================
    # df = trial_log[trial_log.block_id == block_id]
    # all_crit_met = evaluate_block(config,
    #                             df=df,
    #                             block_id = block_id,
    #                             block_log = block_log,
    #                             trial_ons_next_block=trial_ons_next_block)

    return block_log #, all_crit_met



def make_stim_str(stim, framing, m_unit, m_s, m_l, t_unit, t_s, t_l):
    """
    Ths function makes text stimuli to use as choice instruction and choice
    option for the various framings
    :param stim:
    :param m_unit:
    :param m_s:
    :param m_l:
    :param t_unit:
    :param t_s:
    :param t_l:
    :return:
    """

    # Anonymous function for adjusting temporal delay in strings:
    # - replace "in 0 days" with "today"
    # - replace "in 1 days/weeks/months/years" with "in 1
    # day/week/month/year"
    replace_days = \
        lambda choice_str: \
            choice_str. \
                replace("in 0 days", "today"). \
                replace("in 1 " + t_unit, "in 1 " + t_unit[0:-1])

    # We want each line of the choice instructions and choice options to be
    # center-aligned, which we achieve using string formatting:
    if stim == 'choice_instr':
        if framing in ['neutral', 'delay', 'date']:
            stim_str = '{:^70}'.format('Choose between:')
        elif framing in ['defer', 'speedup']:
            fmt_str = "You are entitled to receive {0:s} {1:.2f}  in {2:d} {3:s}. Choose between:"
            if framing == 'defer':
                # Defer framing example:
                # "You are entitled to receive €21.76 today. Choose between:"
                stim_str = \
                    '{:^70}'.format(
                        fmt_str.format(m_unit, m_s, t_s, t_unit))
            elif framing == 'speedup':
                # Speedup framing example:
                # "You are entitled to receive €43.52 in 64 days. Choose between:"
                stim_str = \
                    '{:^70}'.format(
                        fmt_str.format(m_unit, m_l, t_l, t_unit))

        # Adjust string, if needed
        stim_str = replace_days(stim_str)

    elif stim in ['ss_opt', 'll_opt']:

        cal_date_str = \
            lambda t_unit, dt: \
                (calendar.datetime.datetime.today() +
                 {'days': calendar.datetime.timedelta(days=dt),
                  'weeks': calendar.datetime.timedelta(weeks=dt),
                  'months': monthdelta.monthdelta(months=dt),
                  'years': monthdelta.monthdelta(months=12 * dt),
                  }[t_unit]).strftime('%B %d, %Y')

        if stim == 'ss_opt':
            m = m_s
            t = t_s
        elif stim == 'll_opt':
            m = m_l
            t = t_l

        # First row: what to do
        if framing in ['neutral', 'delay', 'date']:
            first_row = '{0:s}'.format('Receive:')
        elif framing == 'defer':
            if stim == 'ss_opt':
                first_row = '{0:s}'.format('As planned, receive:')
            elif stim == 'll_opt':
                first_row = '{0:s}'.format('Defer and receive:')
        elif framing == 'speedup':
            if stim == 'ss_opt':
                first_row = '{0:s}'.format('Speed up and receive:')
            elif stim == 'll_opt':
                first_row = '{0:s}'.format('As planned, receive:')

        # Second row: monetary amount
        second_row = '{0:s} {1:.2f}'.format(m_unit, m)

        # Third row: temporal delay
        if framing == 'date':
            third_row = '{0:s}'.format(cal_date_str(t_unit=t_unit,
                                                    dt=t))
        else:
            third_row = 'in {0:d} {1:s}'.format(t, t_unit)



        third_row = replace_days(third_row)

        # Make the stimulus: a 3-row string center aligned
        stim_str = '{:^70}\n{:^70}\n{:^70}'.format(first_row,
                                                   second_row,
                                                   third_row)

    # Return the stimulus
    return stim_str


def make_trial_list(config):
    """
    Makes a trial list for the indifference point procedure

    :param config:
    :return:
    """

    ip_procedure = config['ip_procedure']

    trial_dur = config['ip_procedure']['trial_dur']
    iti_dur = config['ip_procedure']['iti_dur']

    n_staircase_steps = ip_procedure['n_staircase_trial']
    n_catch_ss = ip_procedure['n_catch_trial_ss']
    n_catch_ll = ip_procedure['n_catch_trial_ll']
    n_instr_check = ip_procedure['n_instr_manip_check_trial']

    n_trial_per_t_l = \
        n_staircase_steps + \
        n_catch_ss + \
        n_catch_ll + \
        n_instr_check
    t_ls = ip_procedure['t_l']

    trial_type = n_catch_ss * ['catch_ss'] + \
                 n_catch_ll * ['catch_ll'] + \
                 n_instr_check * ['instr_check'] + \
                 n_staircase_steps * ['standard']

    trial_list = \
        pd.DataFrame(
            list(it.product(range(n_trial_per_t_l), t_ls)),
            columns = ['trial_ix', 't_l'])
    trial_list['session_ix'] = config['session']['session_ix']
    trial_list['wait_for_trigger'] = False
    trial_list['framing'] = 'neutral'
    trial_list['m_s'] = None
    trial_list['m_l'] = ip_procedure['m_l']
    trial_list['m_units'] = ip_procedure['m_units']
    trial_list['t_s'] = ip_procedure['t_s']
    trial_list['t_units'] = ip_procedure['t_units']
    trial_list['choice_instr_ons'] = 0
    trial_list['choice_instr_dur'] = trial_dur
    trial_list['ss_opt_ons'] = 0
    trial_list['ss_opt_dur'] = trial_dur
    trial_list['ll_opt_ons'] = 0
    trial_list['ll_opt_dur'] = trial_dur

    col_order = ['session_ix', 'block_ix', 'trial_ix', 'trial_ons',
                 'wait_for_trigger', 'framing', 'trial_type', 'm_s', 'm_l',
                 'm_units', 't_s', 't_l', 't_units', 'choice_instr_ons',
                 'choice_instr_dur', 'ss_opt_ons', 'ss_opt_dur',
                 'll_opt_ons', 'll_opt_dur']

    # Each delay is a new block, so participant can take a break between
    # blocks if they wish
    for block_ix, t_l in enumerate(t_ls):
        trial_list.loc[trial_list['t_l'] == t_l, 'block_ix'] = block_ix
        trial_list.loc[trial_list['t_l'] == t_l, 'trial_ons'] = \
            np.arange(0,
                      n_trial_per_t_l * (trial_dur + iti_dur),
                      (trial_dur + iti_dur))
        random.shuffle(trial_type)
        trial_list.loc[trial_list['t_l'] == t_l, 'trial_type'] = trial_type


    # Reorder trial_list's columns and rows
    trial_list = trial_list[col_order]
    trial_list = trial_list. \
        sort_values(by=['t_l', 'trial_ix'], axis='index'). \
        reset_index(drop=True)

    trial_list['trial_ons'] = np.arange(0,
                                        trial_list.shape[0] *
                                        (trial_dur + iti_dur),
                                        (trial_dur + iti_dur))

    return trial_list

def run_ip_procedure(config):
    """
    Runs the indifference point procedure

    Parameters
    ----------
    config      : dict
                specifies itch_time_framing_task experiment properties


    """

    # Retrieve configurations =================================================

    window = config['window']['window']
    stimuli = config['stimuli']
    hub = config['apparatus']['hub']
    rd = config['apparatus']['rd']
    kb = config['apparatus']['kb']
    # trial_stats = config['statistics']['trial']
    ip_procedure = config['ip_procedure']

    m_l = float(ip_procedure['m_l'])
    m_units = ip_procedure['m_units']
    t_l_list = ip_procedure['t_l']
    t_s = int(ip_procedure['t_s'])
    t_units = ip_procedure['t_units']

    m_s = m_l * 0.5


    n_staircase_steps = ip_procedure['n_staircase_trial']
    n_catch_ss = ip_procedure['n_catch_trial_ss']
    n_catch_ll = ip_procedure['n_catch_trial_ll']
    n_instr_check = ip_procedure['n_instr_manip_check_trial']

    # Run Freye et al. indifference point procedure ===========================

    for i_t_l, t_l in enumerate(t_l_list):  # Loop over delays

        t_l = int(t_l)

        m_s = 0

        # for i_trial, trial in enumerate(trials):  # Loop over trials
        for i_trial in range(8):  # Loop over trials

            m_s = m_s + m_l * 2 ** -(i_trial + 1)

            framing = 'neutral'

            #  Specify arguments for run_trial =======

            # Trial list (trial_list) -----------------------------------------
            # TODO: Complete
            trial_dict = {'session_ix': 0,
                          'block_ix': 0,
                          'trial_ix': 0,
                          'trial_ons': 0,
                          'wait_for_trigger': False,
                          'framing': framing,
                          'm_s': 0,
                          'm_l': m_l,
                          'm_units': m_units,
                          't_s': t_s,
                          't_l': t_l,
                          't_units': t_units,
                          'choice_instr_ons': 0,
                          'choice_instr_dur': 10,
                          'ss_opt_ons': 2,
                          'ss_opt_dur': 8,
                          'll_opt_ons': 2,
                          'll_opt_dur': 8
                          # etc.
                          }
            trial_list = pd.DataFrame(trial_dict,index=[0])

            # Trial onset (trial_ons) -----------------------------------------
            trial_ons = 0

            # Trial log (trial_log) -------------------------------------------
            trial_list_ixs = trial_list.index.tolist()
            trial_cols = config['log']['performance']['trial']['columns']
            trial_log = pd.DataFrame(index=trial_list_ixs,
                                     columns=trial_cols)

            # Trial timing (trial_timing) -------------------------------------
            # Mock times, tos that trial starts immediately
            trial_timing = {'ons': -float('inf'),
                            'dur': 0,
                            'iti_dur': 0,
                            'refresh_time': 1 / config['window']['frame_rate']}

            # stim_list, u, f_on_off, trial_dur -------------------------------

            this_trial_log, stim_list, u, f_on_off, trial_dur = \
                stim_to_frame_mat(config,
                                  trial_list.ix[0],
                                  trial_log)

            # Stimuli (stimuli) -----------------------------------------------

            x_ss, y_ss = config['stimuli']['ss_opt'][0].pos
            x_ll, y_ll = config['stimuli']['ll_opt'][0].pos

            # Randomly determine side of SS and LL option
            if random.random() < 0.5:
                stimuli['ss_opt'][0].pos = (x_ss, y_ss)
                stimuli['ll_opt'][0].pos = (x_ll, y_ll)
                this_trial_log.iloc[0]['ll_side'] = 'right'
            else:
                # Swap sides
                stimuli['ss_opt'][0].pos = (x_ll, y_ll)
                stimuli['ll_opt'][0].pos = (x_ss, y_ss)
                this_trial_log.iloc[0]['ll_side'] = 'left'


            for stimulus in ['choice_instr', 'ss_opt', 'll_opt']:
                stimuli[stimulus][0].alignHoriz = 'center'
                stimuli[stimulus][0].alignVert = 'center'
                stimuli[stimulus][0].wrapWidth = 1000

                stimuli[stimulus][0].setText(
                    make_stim_str(stim=stimulus,
                                  framing=framing,
                                  m_unit=m_units,
                                  m_s=m_s,
                                  m_l=m_l,
                                  t_unit=t_units,
                                  t_s=t_s,
                                  t_l=t_l)
                )

            # Run trial
            trial_log = run_trial(config,
                                  wait_for_trigger=False,
                                  trial_ons=trial_ons,
                                  hub=hub,
                                  trial_log=trial_log,
                                  trial_timing=trial_timing,
                                  window=window,
                                  stim_list=stim_list,
                                  u=u,
                                  f_on_off=f_on_off,
                                  rd=rd,
                                  kb=kb,
                                  trial_stats=None,
                                  trial_eval_data=None,
                                  feedback_dur=2,
                                  stimuli=stimuli
                                  )









    m_l = np.double(ip_procedure['m_l'])
    t_ls = np.double(ip_procedure['t_l'])
    t_s = np.double(ip_procedure['t_s'])
    t_units = ip_procedure['t_units']


    m_s = np.around(m_l, decimals = 2) * 0.5


    trials = n_catch_ss * ['catch_ss'] + \
             n_catch_ll * ['catch_ll'] + \
             n_instr_check * ['instr_check'] + \
             n_staircase_steps * ['staircase']

    adjustment_factor = -1


    # Run Freye et al. indifference point procedure ===========================

    for t_l in t_ls: # Loop over delays

        # Shuffle trial sequence, so that catch trials and instruction
        # manipulation check trials, if any, are randomly interspersed
        random.shuffle(trials)

        # Staircase index
        i_staircase = np.cumsum(
            np.int32([trial == 'staircase' for trial in trials]))

        for i_trial, trial in enumerate(trials): # Loop over trials

            # Determine step size
            if trial in ['staircase', 'instr_check']:

                exponent = {'staircase': -np.double(i_staircase[i_trial]),
                            'instr_check': -np.double(i_staircase[i_trial] + 1),
                            }

                step_size = m_l * np.power(np.double(2), exponent[trial])
                m_s = m_s + adjustment_factor * step_size


            if trial == 'catch_ss':
                m_s = m_l
            elif trial == 'catch_ll':
                m_s = 0

            m_l = np.random.rand() * m_l

            # Update stimuli --------------------------------------------------

            # Instruction
            choice_instruction.setText(text='Choose between')
            choice_instruction.setFont(font='Helvetica')
            choice_instruction.setPos(newPos=[0,2])
            choice_instruction.draw()

            # SS option
            ss_opt_str = '%.2f Euro /n in %d days' % (m_s, t_s)
            ss_option = init_stimulus(window, stim_type="TextStim")
            ss_option.setText(text=ss_opt_str )
            choice_instruction.setFont(font='Helvetica')
            ss_option.setPos(newPos=[-8, 0])
            ss_option.draw()

            # SS option
            ll_opt_str = '%.2f Euro /n in %d days' % (m_l, t_l)
            ll_option = init_stimulus(window, stim_type="TextStim")
            ll_option.setText(text=ll_opt_str)
            choice_instruction.setFont(font='Helvetica')
            ll_option.setPos(newPos=[8, 0])
            ll_option.draw()

            window.flip()
            pp.core.wait(0.5)
            window.flip()


            # =================================================================

            trial_ons



            # Run trial
            run_trial(config,
                      wait_for_trigger=True,
                      trial_ons=None,
                      hub=hub,
                      trial_log=None,
                      trial_timing=None,
                      window=window,
                      stim_list=None,
                      u=None,
                      f_on_f_off=None,
                      rd=rd,
                      kb=kb,
                      trial_stats=None,
                      trial_eval_data=None,
                      feedback_dur=None,
                      stimuli=None
                      )





    # def specify_choice_stim(framing, m_ref=0, t_ref=0, dm, dt, m_unit, \
    #                         m_unit_pos='prefix', t_unit, t_unit_pos='suffix'):
    #     """
    #
    #     :param framing:
    #     :param m_ref:
    #     :param t_ref:
    #     :param dm:
    #     :param dt:
    #     :param m_unit:
    #     :param m_unit_pos:
    #     :param t_unit:
    #     :param t_unit_pos:
    #     :return:
    #     """
    #     assert framing in ['neutral', 'delay', 'date', 'defer', 'speedup'], \
    #         'framing should be any of "neutral", "delay", "date", "defer", "speedup"'
    #
    #
    #     bla = {'days': calendar.datetime.timedelta(days=dt),
    #            'weeks': calendar.datetime.timedelta(weeks=dt),
    #            'months': monthdelta.monthdelta(months=dt),
    #            'years': monthdelta.monthdelta(months=12 * dt),
    #            }
    #
    #     t = calendar.datetime.datetime.now()
    #
    #
    #
    #
    #     # Make string
    #     if framing in ['neutral', 'defer', 'speedup', 'delay']:
    #         if m_unit_pos == 'prefix':
    #             choice_str = '%s %.2f in %.0f %s ' % (m_unit,
    #                                                   m_ref + dm,
    #                                                   t_unit,
    #                                                   t_ref + dt)
    #         elif m_unit_pos == 'suffix':
    #             choice_str = '%.2f %s in %.0f %s ' % (m_ref + dm,
    #                                                   m_unit,
    #                                                   t_unit,
    #                                                   t_ref + dt)
    #     elif framing in ['date']:
    #         if m_unit_pos == 'prefix':
    #             choice_str = '%s %.2f in %.0f %s ' % (m_unit,
    #                                                   m_ref + dm,
    #                                                   t_unit,
    #                                                   t_ref + dt)
    #         elif m_unit_pos == 'suffix':
    #             choice_str = '%.2f %s in %.0f %s ' % (m_ref + dm,
    #                                                   m_unit,
    #                                                   t_unit,
    #                                                   t_ref + dt)


    # ss = specify_choice_stim(framing='neutral',
    #                          m_ref=0,
    #                          t_ref=0,
    #                          dm=10,
    #                          dt=0,
    #                          m_unit='Euro',
    #                          m_unit_pos='prefix',
    #                          t_unit='days',
    #                          m_unit_pos='suffix')


    # dt = 3
    # t_units = 'days'
    #
    #
    #
    # my_date = calendar.datetime.datetime.today() +  bla[t_units]
    # my_date.strftime('%B %d, %Y')
    #
    # window = config['window']['window']
    #
    # framing = 'neutral'
    #
    # m_s = 35
    # m_l = 43.52
    # t_s = 0
    # t_l = 4
    # m_units = 'Euro'
    #
    #
    # if framing in ['neutral', 'delay', 'defer', 'speedup']:
    #     fmt_str_choice_opt = '%s %.2f in %.0f %s'
    # elif framing in ['date']:
    #     fmt_str_choice_opt = '%s %.2f in %.0f %s'



    # ss_str = fmt_str_choice_opt % (m_units, m_s, t_s, t_units)
    # ll_str = fmt_str_choice_opt % (m_units, m_s, t_s, t_units)

    # sprintf(    fmt_str_neutral


    # fmt_str_choice_option = {'neutral': "%s %.2f in %.0f %s", m_s, t_s, \
    #                         t_units}


    strFormatPerformance = '%s_Study_%s_TaskVersion_%s_Group_%.2d_Subject_%.3d'


    # m_l = config['ip_procedure']['m_l']
    #
    # task_instruction = pp.visual.TextStim(window,
    #                                       text='Choose between:',
    #                                       pos=(0, 3),
    #                                       color=(1.0, 1.0, 1.0),
    #                                       units='deg',
    #                                       ori=0.0,
    #                                       height=2,
    #                                       antialias=True,
    #                                       bold=False,
    #                                       italic=False,
    #                                       alignHoriz='center',
    #                                       alignVert='center',
    #                                       name='choice_instruction')
    #
    #
    #
    # task_instruction().draw()
    # window.flip()
    # pp.core.wait(2)
    # pp.core.quit()



    # for t_l in config['ip_procedure']['t_l']:


def run_stage(config, stage_id, trial_list):
    """
    Runs all trials comprising the practice, indifference point procedure, or
    experimental stage

    Parameters
    ----------
    config      : dict
                specifies ItchPy experiment properties

    stage_id    : str or unicode
                stage identifier: should be one of 'practice',
                'indifference_point_procedure' or 'experiment'

    trial_list   : pandas.core.frame.DataFrame
                list of trials to present

    """

    block_ixs = trial_list['block_ix'].unique()
    block_cols = config['log']['performance']['block']['columns']
    block_log = pd.DataFrame(index=block_ixs,
                             columns=block_cols)

    sess_columns = config['log']['performance']['sess_columns']
    sess_data = config['log']['performance']['sess_data']
    #
    # performance_req = config['performance_requirements'][stage_id]

    for block_ix in block_ixs:

        all_crit_met = False
        n_iter = 0
        # force_repeat = performance_req['force_repeat']
        # max_n_iter = performance_req['max_n_iter']

        block_id = '%s%.3d' % (stage_id[0], block_ix)

        this_block_log = pd.DataFrame(index=[block_ix], columns=block_cols)

        this_block_log.loc[block_ix, sess_columns] = sess_data
        this_block_log.loc[block_ix, 'block_id'] = block_id
        this_block_log.loc[block_ix, 'block_ix'] = block_ix
        this_block_log.loc[block_ix, 'iterX'] = n_iter

        if block_ix == block_ixs[-1]:
            trial_ons_next_block = np.inf
        else:

            nextblock_ix = [block_ixs[index + 1] for index, value in
                            enumerate(block_ixs) if value == block_ix]
            # TODO: statement seems to have no effect; check if error or
            # what its purpose could be
            # trial_list.loc[trial_list['block_ix'] == nextblock_ix, 'trial_ons']
            trial_ons_next_block = \
                trial_list[trial_list['block_ix'] == nextblock_ix].iloc[0][
                    'trial_ons']

        while not all_crit_met:

            print 'Run %s block %d' % (stage_id, block_ix)

            trial_list_block = trial_list[trial_list['block_ix'] == block_ix]

            # present_instruction(config, stage_id, block_ix)

            # TODO: add output argument allCritMet
            this_block_log = run_block(config=config,
                                       block_id=block_id,
                                       trial_list=trial_list_block,
                                       block_log=this_block_log,
                                       trial_ons_next_block=trial_ons_next_block)

            # Write block log
            # ---------------------------------------------------------
            with open(config['log']['performance']['block']['file'],
                      'a+') as fileObj:
                this_block_log.to_csv(fileObj, index=False, header=False,
                                      na_rep=np.nan)

            if force_repeat:
                if not all_crit_met:
                    if n_iter == (max_n_iter - 1):
                        present_instruction(config, 'block_repeat', 1)
                        present_instruction(config, 'end')
                        pp.core.wait(5)
                        pp.core.quit()
                    else:
                        n_iter = n_iter + 1

                        # Warn subject that block will be repeated
                        present_instruction(config, 'block_repeat', 0)

                        # Reset clock to trial_ons in trial_list
                        # config['clock']['tracking'].reset(trial_list_block.iloc[0]['trial_ons'])
                        config['clock'].reset(
                            trial_list_block.iloc[0]['trial_ons'])
            else:
                break

def run_trial(config,wait_for_trigger,trial_ons,hub,trial_log,trial_timing,window,stim_list,u,f_on_off,rd,kb,trial_stats,trial_eval_data,feedback_dur,stimuli):
    """
    Runs a trial

    Parameters
    ----------
    config          : dict
                    specifies StPy experiment properties

    wait_for_trigger  : bool
                    whether or not trial onsent is contingent on trigger

    trial_ons        : int or float
                    trial onset (in seconds from clock reset)

    hub             : psychopy.iohub.client.ioHubConnection
                    interface to the ioHub

    trial_log        : pandas.core.frame.DataFrame
                    empty trial log

    trial_timing     : dict
                    specifies timing properties, such as trial onset, trial
                    duration, inter-trial interval duration, and refresh times

    window          : psychopy.visual.window.Window
                    PsychoPy window object, in which stimuli are presented

    stim_list       : list
                    list of PsychoPy stimuli

    u               : numpy.ndarray
                    stimulus-by-frame array specifying for each combination of
                    stimulus (row) and frame (column) whether or not the
                    stimulus should be displayed

    f_on_off        : numpy.ndarray
                    stimulus-by-2 array specifying the frame indices of
                    stimulus onset (first column) and stimulus offset (second
                    column)

    rd              : dict
                    specifies response device properties

    kb              : dict
                    specifies keyboard properties

    trial_stats      : dict
                    specifies which descriptive statistics need to be computed

    trial_eval_data   : dict
                    specifies all data necessary to evaluate trials

    feedback_dur     : int or float
                    duration of trial feedback (in seconds)

    stimuli         : dict
                    all stimuli used in the experiment

    Returns
    -------
    trial_log        : pandas.core.frame.DataFrame
                    trial log


    """

    # Wait for external trigger, or trial onset, or both
    # -------------------------------------------------------------------------
    if wait_for_trigger:
        triggered = None
        esc_keys = kb['settings']['esc_keys']
        trigger_keys = kb['settings']['trigger_keys']

        # Clear buffer
        kb['client'].clearEvents()

        # Show text stimulus: Wait for stimulus
        wait_for_trigger_stim = pp.visual.TextStim(window)
        stim_text = u"€" + '  Press any of the following keys to ' \
                    'continue: %s' % \
                    ', '.join(str(x) for x in trigger_keys)
        wait_for_trigger_stim.setText(stim_text)
        wait_for_trigger_stim.draw()
        window.flip()

        # Wait for trigger or abort keys to be pressed
        evs = kb['client'].waitForKeys(keys=esc_keys + trigger_keys)

        if set(esc_keys) & set([ev.key for ev in evs]):
            print('Abort key pressed. Quit PsychoPy now!')
            pp.core.quit()
        elif set(trigger_keys) & set([ev.key for ev in evs]):
            print('Trigger key pressed. Continue now!')
            triggered = True

    if trial_ons == 0:
        # If this is the start of a session or block
        config['clock'].reset()

    while config['clock'].getTime() < (trial_ons - 1.5*trial_timing['refresh_time']):
        pass

    t_start = config['clock'].getTime()

    # Clear events
    # -------------------------------------------------------------------------
    hub.clearEvents('all')

    if __debug__:
        t_events_cleared = config['clock'].getTime()
        print '* Events cleared: t = %f ms; dt = %f ms' % \
              (1000 * (t_events_cleared - t_start),
               1000 * (t_events_cleared - t_start))

    # Present stimuli
    # -------------------------------------------------------------------------
    trial_log = present_stimuli(window=window,
                                stim_list=stim_list,
                                u=u,
                                f_on_off=f_on_off,
                                log=trial_log)

    if __debug__:
        t_stim_presented = config['clock'].getTime()
        print '* Stimuli presented: t = %f ms, dt = %f ms' % \
              (1000*(t_stim_presented-t_start),1000*(t_stim_presented-t_events_cleared))

    # Collect responses
    # -------------------------------------------------------------------------
    trial_log = collect_response(rd=rd, kb=kb, log=trial_log)

    if __debug__:
        t_resp_collected = config['clock'].getTime()
        print '* Responses collected: t = %f ms, dt = %f ms' % \
              (1000*(t_resp_collected-t_start),1000*(t_resp_collected-t_stim_presented))


    # Evaluate trial
    # -------------------------------------------------------------------------
    # trial_log = evaluate_trial(eval_data=trial_eval_data,
    #                            feedback_dur=feedback_dur,
    #                            window=window,
    #                            stimuli=stimuli,
    #                            log=trial_log)
    #
    # if __debug__:
    #     t_feedback_given = config['clock'].getTime()
    #     print '* Feedback given: t = %f ms, dt = %f ms' % \
    #           (1000*(t_feedback_given-t_start),1000*(t_feedback_given-t_stats_computed))

    # Wrap up
    # -------------------------------------------------------------------------
    return trial_log

def stim_to_frame_mat(config,trial,log):
    """
    <SUMMARY LINE>

    <EXTENDED DESCRIPTION>

    Parameters
    ----------
    config      : dict
                specifies StPy experiment properties

    trial       : pandas.core.series.Series
                boolean array acting as selector of trials (rows)

    log         : pandas.core.frame.DataFrame
                trial log

    Returns
    -------
    log         : pandas.core.frame.DataFrame
                trial log; stim_to_frame_mat fills in stimulus properties,
                including information about onsets, durations

    stim_list   : list
                list of PsychoPy stimuli

    u           : numpy.ndarray
                stimulus-by-frame array specifying for each combination of
                stimulus (row) and frame (column) whether or not the stimulus
                should be displayed

    f_on_off    : numpy.ndarray
                stimulus-by-2 array specifying the frame indices of stimulus
                onset (first column) and stimulus offset (second column)

    trial_dur   : numpy.ndarray
                trial duration (in seconds), based on onsets and durations of
                the choice_instr, ss_opt, and ll_opt stimulus.

    """

    trial_log_ix = log.index.tolist()[0]

    stimuli = config['stimuli']

    def append_it(trial,log,stim,stim_list,ons,dur):
        """

        Appends stimulus list

        Parameters
        ----------
        trial       : pandas.core.series.Series
                    trial data

        log         : pandas.core.frame.DataFrame
                    trial log

        stim        : str or unicode
                    stimulus name

        stim_list   : list
                    list of stimuli

        ons         : numpy.ndarray
                    array of stimulus onsets

        dur         : numpy.ndarray
                    array of stimulus durations

        Returns
        -------
        stim_list   : list
                    list of stimuli

        new_ons      : numpy.ndarray
                    array of stimulus onsets

        new_dur      : numpy.ndarray
                    array of stimulus durations

        log         : pandas.core.frame.DataFrame
                    trial log

        """

        # if not pd.isnull(trial[stim + 'Ix']):

        # Identify stimulus index
        # i = trial.loc[stim + 'Ix'].astype(int)

        # Append stimulus list, onset array, and duration array
        stim_list.append(stimuli[stim][0])
        stim_ons = float(trial[stim + '_ons'])
        stim_dur = float(trial[stim + '_dur'])
        new_ons = np.hstack([ons,stim_ons])
        new_dur = np.hstack([dur,stim_dur])

            # log.loc[trial_log_ix,[stim + '_ix']] = i

            # N.B. These are the intended stimulus onsets and durations. They
            # are replace by the actual onsets and durations in present_stimuli
        log.loc[trial_log_ix,[stim + '_ons']] = stim_ons
        log.loc[trial_log_ix,[stim + '_dur']] = stim_dur

        # else:
        #     new_ons = ons
        #     new_dur = dur
        #
        #     for col in ['_ix','_ons','_ons_dt','_dur','_dur_dt']:
        #         log.loc[trial_log_ix,[stim + col]] = np.nan

        return stim_list, new_ons, new_dur, log

    stim_list = []
    ons = np.array([],dtype=int)
    dur = np.array([],dtype=int)

    stim_set = set(config['stimuli'].keys())

    # TODO: Add feedback again
    for stim in stim_set.intersection(['choice_instr', 'ss_opt', 'll_opt']):
        stim_list, ons, dur, log = append_it(trial=trial,
                                             log=log,
                                             stim=stim,
                                             stim_list=stim_list,
                                             ons=ons,
                                             dur=dur)

    # Append iti stimulus to stim_list
    stimuli['iti'][0].setText('o')
    stim_list.append(stimuli['iti'][0])

    dt = np.array([1./config['window']['frame_rate']])
    trial_dur = np.array(np.max(ons + dur))

    # Make stimulus-by-frame matrix (u)
    u, f_on_off, t = time_to_frame(ons=ons, dur=dur, dt=dt, trial_dur=trial_dur)

    return log, stim_list, u, f_on_off, trial_dur

def time_to_frame(ons, dur, dt, trial_dur):
    """
    Makes a stimulus-frame array based on stimulus onsets and durations

    Parameters
    ----------
    ons         : numpy.ndarray
                1D-array of stimulus onset(s)

    dur         : numpy.ndarray
                1D-array of stimulus duration(s)

    dt          : numpy.ndarray
                time step (in seconds)

    trial_dur   : numpy.ndarray
                trial duration (in seconds), based on onsets and durations of
                the choice_instr, ss_opt, and ll_opt stimulus.

    Returns
    -------
    u           : numpy.ndarray
                stimulus-by-frame array specifying for each combination of
                stimulus (row) and frame (column) whether or not the stimulus
                should be displayed

    f_on_off    : numpy.ndarray
                stimulus-by-2 array specifying the frame indices of stimulus
                onset (first column) and stimulus offset (second column)

    t           : numpy.ndarray
                array of frame onset times (in seconds, relative to trial
                onset)

    Raises
    ------
    AssertionError:
        If 'ons' is not instance of class np.ndarray
        If 'dur' is not instance of class np.ndarray
        If 'dt' is not instance of class np.ndarray
        If 'trial_dur' is not instance of class np.ndarray

    """

    ###########################################################################
    # 1. PROCESS INPUTS & SPECIFY VARIABLES
    ###########################################################################

    # 1.1. Import libraries
    # =========================================================================


    # 1.2. Process inputs
    # =========================================================================
    # Check if all inputs are ndarrays

    assert isinstance(ons, np.ndarray), 'ons should be of type ndarray'
    assert isinstance(dur, np.ndarray), 'dur should be of type ndarray'
    assert isinstance(dt, np.ndarray), 'dt should be of type ndarray'
    assert isinstance(trial_dur, np.ndarray), 'trial_dur should be of type ndarray'

    #TODO: assert that inputs are no unsized objects;

    # 1.3. Define dynamic variables
    # =========================================================================
    # Stimulus onsets and offsets, in frames
    f_on = np.around(ons / dt).astype(int)
    f_off = np.around((ons + dur) / dt).astype(int)

    # Number of stimuli
    n_stim = np.size(ons)

    # Number of frames
    if np.isinf(trial_dur):
        n_frame = max(f_off)
    else:
        n_frame = int(np.ceil(trial_dur/dt))

    # Array of frames
    f = range(0, n_frame, 1)

    # Array of time points
    t = np.linspace(0, n_frame * float(dt), n_frame, True)

    # Preallocate u
    u = np.zeros((n_stim, n_frame)).astype(int)

    ###########################################################################
    # 2. SPECIFY STIMULUS AND INPUT MATRICES
    ###########################################################################

    # Loop over stimuli
    for i_stim in range(0, n_stim):

        # Stimulus onsets in stimulus-frame diagram
        u[i_stim, f == f_on[i_stim]] += 1

        # Stimulus offsets in stimulus-frame diagram
        u[i_stim, f == f_off[i_stim]] -= 1

    # Encode 'stimulus on' as 1.0, 'stimulus off' as 0.0
    u = np.cumsum(u, axis=1)

    # Convert to boolean
    u = u.astype(bool)

    # Stimulus on- and offset in frames
    f_on_off = np.vstack([f_on,f_off])

    ########################################
    # 3. SPECIFY OUTPUT
    ########################################
    return u, f_on_off, t