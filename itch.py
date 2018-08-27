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
import warnings

from pprint import pprint
from psychopy.hardware.emulator import launchScan

print 'psychopy version: %s' % ppinfo.psychopyVersion
print 'numpy version: %s' % np.__version__
print 'pandas version: %s' % pd.__version__

# This is to enable printing of all data frame columns
pd.set_option('display.max_columns', None)

def check_df_from_csv_file(df):
    """
    <SUMMARY LINE>

    <EXTENDED DESCRIPTION>

    Parameters
    ----------
    <NAME> : <TYPE>
        <DESCRIPTION>

    Returns
    -------
    <NAME> : <TYPE>
        <DESCRIPTION>

    Raises
    ------
    <EXCEPTIONS>

    Usage
    -----
    <USAGE>

    Example
    -------
    <EXAMPLE THAT CAN IDEALLY BE COPY PASTED>
    """

    # Index (*Ix) and keycount (keycount*) columns should be of type object
    # This is to guarantee that NA and integers can be represented. Floats would
    # cause problems.

    # cols = [col for col in df.select_dtypes(exclude = ['int'])
    #           if col.endswith('Ix') or col.startswith('keyCount')]
    # for col in cols:
    #     df[col] = df[col].astype('object')


    # if os.path.isfile(trialListFile):
    #     trialList = pd.read_csv(trialListFile)
    #     ixCols = [col for col in trialList if re.search('Ix$',col)]
    #
    #
    #     # Assertions
    #     # ---------------------------------------------------------------------
    #     for col in ixCols:
    #         assert trialList[col].dtype == np.int or \
    #                all(trialList['cueIx'].isnull()), \
    #             'column {} in file {} contains data other than integers'.format(col,trialListFile)
    #
    #
    #
    #     config['practice']['trialList'] = trialList

    return df
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

    log         : pandas.core.frame.Series
                trial log

    other_keys  : list (optional)
                specifies which other keys (from response device or keyboard)
                to monitor

    t0          : float (optional)
                event time relative to which RT should be calculated

    min_rt      : float (optional)
                minimum response time, in seconds

    max_rt      : float (optional)
                maximum response time, in seconds

    Returns
    -------
    log         : pandas.core.frame.Series (optional)
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
    other_keys = None
    ev_keys_pressed = None
    other_keys_pressed = None
    log = None

    if len(args) == 0:
        if kwargs:
            if 'log' in kwargs:
                log = kwargs.get('log')
            if 'other_keys' in kwargs:
                other_keys = kwargs.get('other_keys')
            if 't0' in kwargs:
                t0 = kwargs.get('t0')
            else:
                t0 = True
            if 'min_rt' in kwargs:
                min_rt = kwargs.get('min_rt')
            else:
                min_rt = 0
            if 'max_rt' in kwargs:
                max_rt = kwargs.get('max_rt')
            else:
                max_rt = np.inf
    elif len(args) == 1:
        other_keys = args[0]
        t0 = True
        min_rt = 1.5
        max_rt = np.inf
    elif len(args) == 2:
        other_keys = args[0]
        log = args[1]
        t0 = True
        min_rt = 1.5
        max_rt = np.inf
    elif len(args) == 3:
        other_keys = args[0]
        log = args[1]
        t0 = args[2]
        min_rt = 1.5
        max_rt = np.inf
    elif len(args) == 4:
        other_keys = args[0]
        log = args[1]
        t0 = args[2]
        min_rt = args[3]
        max_rt = np.inf
    elif len(args) == 4:
        other_keys = args[0]
        log = args[1]
        t0 = args[2]
        min_rt = args[3]
        max_rt = args[4]
    elif len(args) > 5:
        # TODO: Add error message
        pass

    # Process inputs
    # -------------------------------------------------------------------------
    rsp_keys = []
    for item in rd['settings']['rsp_keys']:
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

    while sum([key_count[key] for key in rsp_keys]) == 0 and \
            pp.core.getTime() < t0 + max_rt:

        rd_events = rd['client'].getEvents()

        for ev in rd_events:

            if rd_class == 'Keyboard':

                # Have any abort keys been pressed?
                if any([re.findall('^'+key+'$', ev.key) for key in esc_keys]):
                    print('Warning: Escape key pressed - Experiment is terminated')
                    pp.core.quit()

                # Have any other keys been pressed (e.g. toggle keys for moving
                # between instruction screens)
                if other_keys:
                    if any([re.findall('^'+key+'$', ev.key) for key in
                            other_keys]):
                        other_keys_pressed = ev.key
                        print other_keys_pressed
                        return key_count, other_keys_pressed

                # If any of the event keys are in event data
                if any([re.findall('^'+key+'$', ev.key) for key in rsp_keys]):
                    key_count[ev.key] += 1
                    if isinstance(log, pd.core.series.Series):
                        # Only log time of first response key event
                        if not any([key_time[key] for key in rsp_keys]):
                            rt = ev.time - t0
                            key_time[ev.key] = rt
                            break

            elif rd_class == 'Serial':

                for kbev in kb['client'].getEvents():

                    # Have any abort keys been pressed?
                    if any([re.findall('^'+key+'$', kbev.key) for key in esc_keys]):
                        print(
                            'Warning: Escape key pressed - Experiment is '
                            'terminated')
                        pp.core.quit()

                    if other_keys:
                        if any([re.findall('^'+key+'$', kbev.key) for key in
                                other_keys]):
                            other_keys_pressed = kbev.key
                            return key_count, other_keys_pressed

                    # If any of the event keys are in event data
                    if any([re.findall('^'+key+'$', ev.data) for key in rsp_keys]):
                        key_count[ev.data] += 1

                        # Only log time of first response key event
                        if not any([key_time[key] for key in rsp_keys]):
                            rt = ev.time - t0
                            key_time[ev.key] = rt
                            break

    # For each key, only response times of first two events are stored
    if isinstance(log, pd.core.series.Series):


        # # Determine choice
        # try:
        #     choice_key, choice_time = \
        #         next((k, v) for k, v in key_time.items() if v)
        # except StopIteration:
        #     print 'The generator was empty'
        #     choice_key = None
        #     choice_time = float('inf')

        if log['ll_side'] == 'left':
            key_mapping = dict(zip(rsp_keys, ('ll', 'ss', 'avoid_choice')))
        elif log['ll_side'] == 'right':
            key_mapping = dict(zip(rsp_keys, ('ss', 'll', 'avoid_choice')))

        # Log events
        for key in rsp_keys:
            log['key_count_' + key] = key_count[key]
        if 'ev' in locals():
            log['choice'] = key_mapping.get(ev.key, 'NA')
            log['rt'] = rt
            if rt < min_rt:
                log['too_fast'] = True
            else:
                log['too_fast'] = False
        else:
            # If no response is given, variable ev does not exist
            log['choice'] = 'NA'
            log['rt'] = np.inf
            log['too_fast'] = False

        return log

    else:

        return key_count, other_keys_pressed
def copy_series_values(var_names, source_series, target_series):
    assert isinstance(source_series, pd.core.frame.Series), \
        'source_series should be of type pandas.core.frame.Series'
    assert isinstance(target_series, pd.core.frame.Series), \
        'target_series should be of type pandas.core.frame.Series'

    for item in var_names:
        target_series[item] = source_series[item]

    return target_series
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

    if 'content' in stim_info:
        n_stimulus = len(stim_info['content'])
    else:
        n_stimulus = 1

        # Initialize list of stimuli
    stimulus = [None] * n_stimulus

    for i in range(n_stimulus):

        stimulus[i] = init_stimulus(window,stim_info['type'])

        assert (type(stimulus[i]) is
                pp.visual.text.TextStim or
                pp.visual.rect.Rect or
                pp.visual.image.ImageStim), \
                "stimulus is neither a TextStim, nor a Rect, nor an ImageStim"

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


            # Set other parameters, if provided
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

            # Set general TextStimulus parameters
            stimulus[i].alignHoriz = 'center'
            stimulus[i].alignVert = 'center'
            stimulus[i].wrapWidth = 1000

        # Rectangles
        # ---------------------------------------------------------------------
        elif type(stimulus[i]) is pp.visual.Rect:


            if 'height' in stim_info:
                if stim_info['height']:
                    stimulus[i].setHeight(stim_info['height'])

            if 'width' in stim_info:
                if stim_info['width']:
                    stimulus[i].setWidth(stim_info['width'])
            if 'pos' in stim_info:
                if stim_info['pos']:
                    stimulus[i].setPos(stim_info['pos'])
            if 'ori' in stim_info:
                if stim_info['ori']:
                    stimulus[i].setOri(stim_info['ori'])

            if 'fill_color' in stim_info:
                if stim_info['fill_color']:
                    stimulus[i].setFillColor(stim_info['fill_color'], 'rgb255')

            if 'line_color' in stim_info:
                if stim_info['line_color']:
                    stimulus[i].setLineColor(stim_info['line_color'], 'rgb255')

            if 'line_width' in stim_info:
                if stim_info['line_width']:
                    stimulus[i].setLineWidth(stim_info['line_width'])

            if 'opacity' in stim_info:
                if stim_info['opacity']:
                    stimulus[i].setOpacity(stim_info['opacity'])



        # Image stimuli
        # ---------------------------------------------------------------------
        elif type(stimulus[i]) is pp.visual.ImageStim:

            # Set stimulus content
            stimulus[i].setImage(stim_info['content'][i])

            # Set other parameters
            if 'ori' in stim_info:
                if stim_info['ori']:
                    stimulus[i].setOri(stim_info['ori'])

            if 'pos' in stim_info:
                if stim_info['pos']:
                    stimulus[i].setPos(stim_info['pos'])


    return stimulus
def evaluate_block(config,df,block_id,block_log):
    """
    Evaluate block performance

    Parameters
    ----------
    config              : dict
                        specifies ItchPy experiment properties

    df                  : pandas.core.frame.DataFrame
                        trial log

    block_id             : str or unicode
                        identifier of the block

    block_log            : pandas.core.frame.DataFrame
                        block log

    Returns
    -------
    all_crit_met          : bool
                        whether or not all predefined task performance
                        criteria
                        have been met
    """

    # Subfunctions
    # =========================================================================
    def assess_performance(stat, *args, **kwargs):
        """
        Evaluates whether statistic is within lower and upper bounds

        Parameters
        ----------
        stat    : dict
                statistical value to be evaluated

        lo      : int or float
                lower bound

        hi      : int or float
                upper bound

        Returns
        -------
        crit_met : bool
                specifies whether or not criterion is met

        """

        if len(args) == 0:
            if kwargs:
                if 'lo' in kwargs:
                    lo = kwargs.get('lo')
                if 'hi' in kwargs:
                    hi= kwargs.get('hi')
                if 'max_dev' in kwargs:
                    max_dev = kwargs.get('max_dev')
        elif len(args) == 1:
            max_dev = args[0]
        elif len(args) == 2:
            lo = args[0]
            hi = args[1]

        if 'max_dev' in kwargs:
            if all(stat <= max_dev):
                crit_met = True
            else:
                crit_met = False
        else:
            if lo <= stat <= hi:
                crit_met = True
            else:
                crit_met = False

        return crit_met
    def get_bounds(config, stat):
        """
        Get lower and upper bounds for a given statistic

        Parameters
        ----------
        config  : dict
                specifies StPy experiment properties

        stat    : str or unicode
                descriptive statistic name

        Returns
        -------
        lo      : int or float
                lower bound of criterion

        hi      : int or float
                upper bound of criterion

        """

        crit = config['feedback']['block']['features'][stat]['criterion']

        return min(crit), max(crit)
    def get_data(df, stat_type):
        """
        Obtain data from DataFrame based on which descriptive statistics
        are computed

        Parameters
        ----------
        df          : pandas.core.series.DataFrame
                    trial log

        stat_type    : str or unicode
                    statistic on which feedback is presented

        Returns
        -------
        data        : pandas.core.series.Series
                    data on which descriptive statistic are computed

        """

        tt = {'catch_ss_accuracy': 'catch_ss',
              'catch_ll_accuracy': 'catch_ll',
              'check_instr_accuracy': 'check_instr',
              'too_fast_responses': 'standard',
              'monotonicity': 'standard',
              'discounting': 'standard'
              }[stat_type]


        if stat_type.endswith('accuracy'):
            data = df.trial_correct[df.trial_type == tt].value_counts()
        elif stat_type in ['monotonicity', 'discounting']:
            # Indifference points
            data = \
                df.loc[df.trial_type == tt, ['trial_ix','t_l', 'm_s']]. \
                sort(columns=['trial_ix', 't_l']). \
                groupby(by='t_l')['m_s']. \
                agg('last')
        elif stat_type == 'too_fast_responses':
            data = df.too_fast[df.trial_type == tt].value_counts()
        else:
            data = None

        return data
    def get_desc_stat(config, data, stat_type):
        """
        Compute descriptive statistic on data

        Parameters
        ----------
        config  : dict
                specifies ItchPy experiment properties

        data        : int or float
                    data on which descriptive statistic is computed

        stat_type    : str or unicode
                    statistic on which feedback is presented

        Returns
        -------
        desc_stat    : int or float
                    descriptive statistic

        """

        # Assertions
        known_stat_types = ['catch_ss_accuracy', 'catch_ll_accuracy',
                            'check_instr_accuracy', 'too_fast_responses',
                            'monotonicity', 'discounting']
        assert stat_type in known_stat_types, 'unknown stat_type %s' % \
                                              stat_type

        if stat_type.endswith('accuracy') or stat_type == 'too_fast_responses':

            if True in data.index:
                n_true = data[True].astype(float)
                n_trial = data.sum().astype(float)
                p_correct = n_true / n_trial
                desc_stat = (p_correct * 100).round()
            else:
                desc_stat = 0

        elif stat_type == 'monotonicity':
            if data.empty:
                desc_stat = np.nan
            else:
                # The differences in indifference point across delays are
                # computed; an indifference point cannot be higher than the
                # previous indifference point by 20% of the LL amount.
                # This is criterion 1 from Johnson & Bickel (2008)
                desc_stat = np.diff(data) / config['ip_procedure']['m_l'] * 100
        elif stat_type == 'discounting':
            if data.empty:
                desc_stat = np.nan
            else:
                # The indifference point of the longest delay is used to
                # evaluate whether indifference points decrease with delay
                # at all
                # This is criterion 2 from Johnson & Bickel (2008)
                desc_stat = data.iloc[-1] / config['ip_procedure']['m_l'] * 100
        else:
            desc_stat = np.nan

        return desc_stat
    def get_feedback_message(config, stat):
        """
        Determine what feedback should be presented

        Parameters
        ----------
        config      : dict
                    specifies StPy experiment properties

        stat        : str or unicode
                    descriptive statistic name

        Returns
        -------
        pos_mes      : str, unicode, or list
                    feedback message if performance criterion is met

        neg_mes      : str, unicode, or list
                    feedback message if performance criterion is not met

        """

        pos_mes = config['feedback']['block']['features'][stat][
            'feedbackPos']
        neg_mes = config['feedback']['block']['features'][stat][
            'feedbackNeg']

        pos = pos_mes
        neg = neg_mes

        return str(pos), str(neg)
    def update_feedback_log(log, stat, stat_type, crit_met):
        """
        Update block log

        Logs performance level and whether or not preset performance
        criterion is met.

        Parameters
        ----------
        log         : pandas.core.frame.DataFrame
                    block log

        stat        : numpy.float64 or numpy.int
                    performance

        stat_type    : str or unicode
                    statistic on which feedback is presented

        crit_met     : bool
                    whether or not performance criterion is met

        Returns
        -------
        log         : pandas.core.frame.DataFrame
                    block log

        """

        # Dict of formatted strings, referring to columns in log
        str_stat_col = {'catch_ss_accuracy': 'catch_ss_accuracy',
                      'catch_ll_accuracy': 'catch_ll_accuracy',
                      'check_instr_accuracy': 'check_instr_accuracy',
                      'too_fast_responses': 'too_fast_responses',
                      'monotonicity': 'monotonicity',
                      'discounting': 'discounting'
                      }
        str_crit_col = {'catch_ss_accuracy': 'catch_ss_accuracy_crit_met',
                      'catch_ll_accuracy': 'catch_ll_accuracy_crit_met',
                      'too_fast_responses': 'too_fast_responses_crit_met',
                      'monotonicity': 'monotonicity_crit_met',
                      'discounting': 'discounting_crit_met'
                      }

        # Column names for statistic and criterion
        col_stat = str_stat_col[stat_type]
        col_crit = str_crit_col[stat_type]

        # Update the log
        if col_stat == 'monotonicity':
            log[col_stat] = max(stat)
        else:
            log[col_stat] = stat
        log[col_crit] = crit_met

        return log
    def update_feedback_screen(win, feedback_stim, stim, stat, stat_type,
                               crit_met, pos_mes, neg_mes):
        """
        Update feedback stimulus

        Parameters
        ----------
        win             : psychopy.visual.window.Window
                        PsychoPy window object, in which stimuli are
                        presented

        feedback_stim    : dict
                        Specifies aspects of the feedback: stimulus
                        identity,
                        performance, and feedback message

        stim            : psychopy.visual.text.TextStim or
        psychopy.visual.text.ImageStim
                        PsychoPy stimulus to which feedback relates

        stat            : numpy.float64 or numpy.int
                        performance

        stat_type        : str or unicode
                        statistic on which feedback is presented

        crit_met         : bool
                        whether or not performance criterion is met

        pos_mes          : str or unicode
                        feedback message if performance criterion is met

        neg_mes          : str or unicode
                        feedback message if performance criterion is not
                        met

        Returns
        -------
        feedback_stim    : dict
                        Specifies aspects of the feedback: stimulus
                        identity,
                        performance, and feedback message

        """

        # Define some variables
        # -----------------------------------------------------------------
        # stimName = stim.name[stim.name.find('_') + 1:]

        stim_name_text = get_empty_text_stim(win)
        perform_text = get_empty_text_stim(win)
        feedback_text = get_empty_text_stim(win)

        pos_feedback_color = (0, 191, 0)
        neg_feedback_color = (255, 128, 0)

        # Stimulus
        # -----------------------------------------------------------------
        stim_str = {'catch_ss_accuracy': "accuracy on catch scenarios with a €0.00 option",
                    'catch_ll_accuracy': 'accuracy on catch scenarios with equal amount options',
                    'check_instr_accuracy': 'accuracy on catch scenarios requiring return key press',
                    'too_fast_responses': 'proportion of fast responses',
                    'monotonicity': 'max. difference between indifference '
                                    'points (%)',
                    'discounting': 'degree of discounting at longest delay'
                    }


        stim_name_text.setText(stim_str[stat_type])
        feedback_stim['stim'].append(stim_name_text)

        # Performance
        # -----------------------------------------------------------------

        stat_str = {'catch_ss_accuracy': '%0.f%%',
                    'catch_ll_accuracy': '%0.f%%',
                    'check_instr_accuracy': '%0.f%%',
                    'too_fast_responses': '%0.f%%',
                    'monotonicity': '%0.f%%',
                    'discounting': '%0.f%%'
                    }

        if stat_type == 'monotonicity':
            perform_text.setText(stat_str[stat_type] % max(stat))
        else:
            perform_text.setText(stat_str[stat_type] % stat)

        if crit_met:
            perform_text.setColor(pos_feedback_color, 'rgb255')
        else:
            perform_text.setColor(neg_feedback_color, 'rgb255')

        feedback_stim['performance'].append(perform_text)

        # Feedback
        # -----------------------------------------------------------------
        if crit_met:
            feedback_text.setText(pos_mes)
            feedback_text.setColor(pos_feedback_color, 'rgb255')
        else:
            feedback_text.setText(neg_mes)
            feedback_text.setColor(neg_feedback_color, 'rgb255')

        feedback_stim['feedback'].append(feedback_text)

        return feedback_stim

    # Criteria
    # =========================================================================

    # Practice
    # -------------------------------------------------------------------------
    # - Responses on all trials

    # Indifference point procedure
    # -------------------------------------------------------------------------
    # - Responses on all trials
    # - Accuracy on catch trials

    # trial_types = ['catch_ss', 'catch_ll', 'instr_check']
    # accuracies = {'catch_ss': None,
    #               'catch_ll': None,
    #               'instr_check': None}




    # #########################################################################

    window = config['window']['window']
    trial_stats = config['statistics']['trial']

    tt = df.trial_type
    trial_type_exist = {'catch_ss_accuracy': any(tt == 'catch_ss'),
                        'catch_ll_accuracy': any(tt == 'catch_ll'),
                        'check_instr_accuracy': any(tt == 'check_instr'),
                        'too_fast_responses': any(tt == 'standard'),
                        'monotonicity': any(tt == 'standard'),
                        'discounting': any(tt == 'standard')
                        }

    # Stimulus
    stimulus = {'catch_ss_accuracy': 'catch_ss',
                'catch_ll_accuracy': 'catch_ll',
                'check_instr_accuracy': 'check_instr',
                'too_fast_responses': 'too_fast_responses',
                'monotonicity': 'ip_decrease',
                'discounting': 'discounting'}

    # Task performance features to provide feedback on
    features = config['feedback']['block']['features']
    feedback_feat = [key for key in features.keys() if
                    features[key]['enable'] and trial_type_exist[key]]

    block_feedback = config['feedback']['block']['features']
    block_feedback_stim = {'stim': [],
                           'performance': [],
                           'feedback': []}

    criteria_met = []

    for feat in sorted(feedback_feat):

        if feat == 'monotonicity':
            max_dev = config['feedback']['block']['features'][feat]['criterion']
        else:
            lower_bound, upper_bound = get_bounds(config=config,
                                                  stat=feat)

        pos_message, neg_message = get_feedback_message(config=config,
                                                        stat=feat)

        data = get_data(df=df,
                        stat_type=feat)

        if not data.empty:
            desc_stat = get_desc_stat(config=config,
                                      stat_type=feat,
                                      data=data)
            if feat == 'monotonicity':
                this_crit_met = assess_performance(stat=desc_stat,
                                                   max_dev=max_dev)
            else:
                this_crit_met = assess_performance(stat=desc_stat,
                                                   lo=lower_bound,
                                                   hi=upper_bound)

            criteria_met.append(this_crit_met)

            # Update feedback screen
            # TODO: implement this
            block_feedback_stim = update_feedback_screen(win=window,
                                                         feedback_stim=block_feedback_stim,
                                                         stim=stimulus[feat],
                                                         stat=desc_stat,
                                                         stat_type=feat,
                                                         crit_met=this_crit_met,
                                                         pos_mes=pos_message,
                                                         neg_mes=neg_message)

            # Update feedback log
            block_log = update_feedback_log(log=block_log,
                                            stat=desc_stat,
                                            stat_type=feat,
                                            crit_met=this_crit_met)

    all_crit_met = all(criteria_met)

    # Display feedback
    # -------------------------------------------------------------------------

    # Count how lines feedback
    n_lines = len(block_feedback_stim['stim'])

    # Feedback title, containing block ID
    block_title_stim = get_empty_text_stim(window)
    y_pos = (float(n_lines) - 1) / 2 + 2
    x_pos = 0
    block_title_stim.setText('Block %s' % (block_id))
    block_title_stim.setPos((x_pos, y_pos))
    block_title_stim.setHeight(1)
    block_title_stim.alignHoriz = 'center'
    block_title_stim.setAutoDraw(True)

    # Loop over feedback lines
    for i_stim in range(n_lines):
        # Set position of the stimulus
        y_pos = (float(n_lines) - 1) / 2 - i_stim
        x_pos = -22.5

        block_feedback_stim['stim'][i_stim].setPos((x_pos, y_pos))
        block_feedback_stim['stim'][i_stim].setHeight(0.75)
        block_feedback_stim['stim'][i_stim].alignHoriz = 'left'
        block_feedback_stim['stim'][i_stim].setAutoDraw(True)

        # Set position of performance stimulus
        x_pos = -2.5
        block_feedback_stim['performance'][i_stim].setPos((x_pos, y_pos))
        block_feedback_stim['performance'][i_stim].setHeight(0.75)
        block_feedback_stim['performance'][i_stim].alignHoriz = 'left'

        block_feedback_stim['performance'][i_stim].setAutoDraw(True)

        # Set position of feedback stimulus
        x_pos = 5
        block_feedback_stim['feedback'][i_stim].setPos((x_pos, y_pos))
        block_feedback_stim['feedback'][i_stim].setHeight(0.75)
        block_feedback_stim['feedback'][i_stim].alignHoriz = 'left'

        block_feedback_stim['feedback'][i_stim].setAutoDraw(True)

    window.flip()

    t_now = pp.core.getTime()
    feedback_duration = config['feedback']['block']['duration']
    pp.core.wait(feedback_duration)

    block_title_stim.setAutoDraw(False)
    for i_stim in range(n_lines):
        block_feedback_stim['stim'][i_stim].setAutoDraw(False)
        block_feedback_stim['performance'][i_stim].setAutoDraw(False)
        block_feedback_stim['feedback'][i_stim].setAutoDraw(False)
    window.flip()

    return all_crit_met
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

                response_type    : pandas.core.series.Series
                                response type for each of the stimulus-response
                                combinations

                feedback        : pandas.core.series.Series
                                feedback for each of the stimulus-response
                                combinations

                trial_type       : pandas.core.series.Series
                                trial type for each of the stimulus-response
                                combinations

    feedback_dur : float
                trial feedback duration (in seconds)

    window      : psychopy.visual.window.Window
                PsychoPy window object, in which stimuli are presented

    stimuli     : dict
                specifies PsychoPy stimuli, including the feedback and inter-
                trial interval stimulus

    log         : pandas.core.frame.Series
                trial log

    Returns
    -------
    log         : pandas.core.frame.Series
                trial log; evaluate_trial fills in values for the following
                variables: trialCorrect, trial_type, response_type, and
                trial_feedback

    """

    # Process inputs
    # =========================================================================

    # Assertions
    # -------------------------------------------------------------------------

    # specific to eval_data and log (window and stimulis should have been
    # checked already and not have been updated

    # Define dynamic variables
    # -------------------------------------------------------------------------
    trial_feedback = []


    trial_correct    = []
    trial_feedback   = []

    # Trial evaluation
    # =========================================================================

    # Match pattern, using stimulus and response data
    # TODO: Only fillna for source_columns, not for e.g. feedback
    source = eval_data['eval_data'].fillna('NA')
    pattern_dict= log[source.columns].fillna('NA').to_dict()
    pattern = {key: [value] for key, value in pattern_dict.iteritems()}
    i_row = source.isin(pattern).all(1)

    if sum(i_row) == 0:
        warnings.warn('No match between trial performance and trial '
                      'evaluation criteria', UserWarning)
        trial_correct = False
        trial_feedback = ""

    elif sum(i_row) == 1:

        trial_correct = \
            eval_data['correct'].get(eval_data['correct'][i_row].index[0])
        trial_feedback = \
            eval_data['feedback'].get(eval_data['feedback'][i_row].index[0])

    log['trial_correct'] = trial_correct
    log['trial_feedback'] = trial_feedback

    # Staircase adjustments
    # =========================================================================

    # Feedback
    # =========================================================================

    stimuli['feedback'][0].setText(trial_feedback)

    if 'll' in log['choice']:
        feedback_stim_to_display = \
            ['feedback', 'll_rect', 'choice_instr', 'ss_opt', 'll_opt']
    elif 'ss' in log['choice']:
        feedback_stim_to_display = \
            ['feedback', 'ss_rect', 'choice_instr', 'ss_opt', 'll_opt']
    else:
        feedback_stim_to_display = \
            ['feedback', 'choice_instr', 'ss_opt', 'll_opt']

    # Display feedback stimuli
    for stim in feedback_stim_to_display :
        stimuli[stim][0].draw()

    window.flip()
    pp.core.wait(feedback_dur)

    return log
def get_empty_text_stim(window):
    """
    Returns an empty PsychoPy text stimulus object

    This object has the following properties:

    text        = '',
    font        = 'Arial',
    pos         = (0,0),
    color       = [1.0, 1.0, 1.0],
    colorSpace  = 'rgb',
    opacity     = 1.0,
    bold        = False,
    italic      = False,
    alignHoriz  = 'center',
    alignVert   = 'center',
    wrapWidth   = 1000,
    autoLog     = None


    Parameters
    ----------
    window      : psychopy.visual.window.Window
                PsychoPy window object, in which stimuli are presented

    Returns
    -------
    textStim    : psychopy.visual.TextStim
                PsychoPy text stimulus object

    """
    text_stim = pp.visual.TextStim(window,
                                   text='',
                                   font='Arial',
                                   pos=(0,0),
                                   color=[1.0, 1.0, 1.0],
                                   colorSpace='rgb',
                                   opacity=1.0,
                                   bold=False,
                                   italic=False,
                                   alignHoriz='center',
                                   alignVert='center',
                                   wrapWidth=1000,
                                   autoLog=None)

    text_stim.setSize(2,units='deg')

    return text_stim
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

    trial_eval_data = pd.read_csv(
        config['evaluation']['trial']['eval_data_file'][rd_class],
        error_bad_lines=False, keep_default_na=False)
    # trial_eval_data = check_df_from_csv_file(trial_eval_data)

    trial_categories = trial_eval_data.fillna(value="")

    eval_columns = [col for col in trial_eval_data.columns if col not in
                  ['correct', 'feedback']]

    config['evaluation']['trial']['eval_data'] = trial_eval_data[
        eval_columns].copy()
    config['evaluation']['trial']['trial_type'] = trial_eval_data[
        'trial_type'].copy()
    config['evaluation']['trial']['choice'] = trial_eval_data[
        'choice'].copy()
    config['evaluation']['trial']['correct'] = trial_eval_data[
        'correct'].copy()
    config['evaluation']['trial']['feedback'] = trial_eval_data[
        'feedback'].copy()

    ###########################################################################
    # INSTRUCTION
    ###########################################################################

    if user_def_params['instruction']:
        config['instruction']['enable'] = True
    else:
        config['instruction']['enable'] = False


    config['stimuli']['instruction'] = \
        {instr_type: define_stimulus(window,
                                     config['instruction'][instr_type],
                                     config['instruction']['experiment'][
                                         'name']
                                     )
         for instr_type in ['start', 'practice', 'ip_procedure',
                       'block_repeat', 'experiment', 'break', 'end']
         if config['instruction'][instr_type]['content'] is not None}

    ###########################################################################
    # PRACTICE
    ###########################################################################

    if user_def_params['practice']:
        config['practice']['enable'] = True
    else:
        config['practice']['enable'] = False


    ###########################################################################
    # INDIFFERENCE POINT PROCEDURE
    ###########################################################################

    if user_def_params['ip_procedure']:
        config['ip_procedure']['enable'] = True
    else:
        config['ip_procedure']['enable'] = False

    config['ip_procedure']['m_l'] = \
        2 ** (config['ip_procedure']['n_staircase_trial'] + 1) * \
        config['ip_procedure']['step_size']

    ###########################################################################
    # EXPERIMENT
    ###########################################################################

    if user_def_params['experiment']:
        config['experiment']['enable'] = True
    else:
        config['experiment']['enable'] = False

    ###########################################################################
    # PERFORMANCE REQUIREMENTS
    ###########################################################################

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
                      'iter_ix',
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
                        'trial_dur']

        # Trigger
        # =====================================================================
        trigger_columns = ['waited_for_trigger']

        # Feedback
        # =====================================================================
        # TODO: Think about feedback logging
        feedback_columns = ['trial_correct',
                            'too_fast',
                            'trial_feedback']
        # feedback_columns = ['trial_correct',
        #                     'trial_type',
        #                     'response_type',
        #                     'trial_feedback']

        # Choice options
        # =====================================================================
        choice_opt_columns = ['trial_type', 'm_s_cat', 'm_s', 'm_l', 'm_unit',
                              't_s', 't_l', 't_unit']

        # Response event times
        # =====================================================================
        rsp_keys = [item
                    for item
                    in config['apparatus']['rd']['settings']['rsp_keys']]

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
        columns += feedback_columns

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
                      'iter_ix']

        # Accuracy - catch-ss
        # =====================================================================
        if config['feedback']['block']['features']['catch_ss_accuracy'][
            'enable']:
            catch_ss_accuracy_cols = ['catch_ss_accuracy',
                                      'catch_ss_accuracy_crit_met']
        else:
            catch_ss_accuracy_cols = []

        # Accuracy - catch-ll
        # =====================================================================
        if config['feedback']['block']['features']['catch_ll_accuracy'][
            'enable']:
            catch_ll_accuracy_cols = ['catch_ll_accuracy',
                                      'catch_ll_accuracy_crit_met']
        else:
            catch_ll_accuracy_cols = []

        # Accuracy - check_instr
        # =====================================================================
        if config['feedback']['block']['features']['check_instr_accuracy'][
            'enable']:
            check_instr_accuracy_cols = ['check_instr_accuracy',
                                         'check_instr_accuracy_crit_met']
        else:
            check_instr_accuracy_cols = []

        # Indifference points increasing monotonically
        # =====================================================================
        if config['feedback']['block']['features']['monotonicity'][
            'enable']:
            monotonicity_cols = ['monotonicity',
                                          'monotonicity_crit_met']
        else:
            monotonicity_cols = []

        # Fast responses
        # =====================================================================
        if config['feedback']['block']['features']['too_fast_responses'][
            'enable']:
            too_fast_responses_cols = ['too_fast_responses',
                                       'too_fast_responses_crit_met']
        else:
            too_fast_responses_cols = []

        # Indifference points decreasing sufficiently
        # =====================================================================
        if config['feedback']['block']['features']['discounting'][
            'enable']:
            discounting_cols = ['discounting',
                                   'discounting_crit_met']
        else:
            discounting_cols = []

        # Put it together
        # =====================================================================
        columns += id_columns
        columns += ix_columns
        columns += catch_ss_accuracy_cols
        columns += catch_ll_accuracy_cols
        columns += check_instr_accuracy_cols
        columns += too_fast_responses_cols
        columns += monotonicity_cols
        columns += discounting_cols

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
                 'rect': pp.visual.Rect(window),
                 'imagestim': pp.visual.ImageStim(window),
                 }

    stim_object = stim_dict[stim_type.lower()]

    return stim_object
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

    assert any([isinstance(stim, str), isinstance(stim, unicode)]), \
        'stim should be of type srr or unicode'
    assert any([isinstance(framing, str), isinstance(framing, unicode)]), \
        'framing should be of type srr or unicode'
    assert any([isinstance(m_unit, str), isinstance(m_unit, unicode)]), \
        'm_unit should be of type srr or unicode'
    assert any([isinstance(m_s, float), isinstance(m_s, int)]), \
        'm_s should be of type float or int'
    assert any([isinstance(m_l, float), isinstance(m_l, int)]), \
        'm_l should be of type float or int'
    assert any([isinstance(t_unit, str), isinstance(t_unit, unicode)]), \
        't_unit should be of type srr or unicode'
    assert any([isinstance(t_s, float), isinstance(t_s, int)]), \
        't_s should be of type float or int'
    assert any([isinstance(t_l, float), isinstance(t_l, int)]), \
        't_l should be of type float or int'

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
            fmt_str = "You are scheduled to receive {0:s} {1:.2f}  in {2:d} {3:s}. Choose between:"
            if framing == 'defer':
                # Defer framing example:
                # "You are scheduled to receive €21.76 today. Choose between:"
                stim_str = \
                    '{:^70}'.format(
                        fmt_str.format(m_unit, m_s, t_s, t_unit))
            elif framing == 'speedup':
                # Speedup framing example:
                # "You are scheduled to receive €43.52 in 64 days. Choose between:"
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
            # first_row = '{0:s}'.format('Receive:')
            first_row = '{0:s}'.format('')
        elif framing == 'defer':
            if stim == 'ss_opt':
                # first_row = '{0:s}'.format('As planned, receive:')
                first_row = '{0:s}'.format('As planned:')
            elif stim == 'll_opt':
                # first_row = '{0:s}'.format('Defer and receive:')
                first_row = '{0:s}'.format('Defer:')
        elif framing == 'speedup':
            if stim == 'ss_opt':
                # first_row = '{0:s}'.format('Speed up and receive:')
                first_row = '{0:s}'.format('Speed up:')
            elif stim == 'll_opt':
                # first_row = '{0:s}'.format('As planned, receive:')
                first_row = '{0:s}'.format('As planned:')

        # Second row: monetary amount
        try:
            second_row = '{0:s} {1:.2f}'.format(m_unit, m)
        except ValueError:
            print('Oops! Something went wrong!')

        # Third row: temporal delay
        if framing == 'date':
            third_row = '{0:s}'.format(cal_date_str(t_unit=t_unit,
                                                    dt=t))
        else:
            third_row = 'in {0:d} {1:s}'.format(t, t_unit)

        third_row = replace_days(third_row)

        # Make the stimulus: a 3-row string center aligned
        if framing in ['neutral', 'delay', 'date']:
            stim_str = '{:^70}\n{:^70}'.format(second_row,
                                               third_row)
        else:
            stim_str = '{:^70}\n{:^70}\n{:^70}'.format(first_row,
                                                       second_row,
                                                       third_row)

    # Return the stimulus
    return stim_str
def make_trial_list(config, stage_id):
    """
    Makes a trial list for the indifference point procedure

    :param config:
    :param stage_id:
    :return:
    """

    # Get variables from config file
    # =========================================================================

    # Variables that are common between stages
    # -------------------------------------------------------------------------

    # Variables related to timing
    self_paced = config[stage_id]['self_paced']
    isi_dur = config[stage_id]['isi_dur']
    max_rt = config[stage_id]['max_rt']
    min_feedback_dur = config['feedback']['trial']['min_duration']
    iti_dur = config[stage_id]['iti_dur']


    if stage_id == 'experiment':
        time_frames = config[stage_id]['time_frames']

    # Variables related to catch trials
    if stage_id in ['ip_procedure', 'experiment']:
        n_catch_ss = config[stage_id]['n_catch_trial_ss']
        n_catch_ll = config[stage_id]['n_catch_trial_ll']
        n_instr_check = config[stage_id]['n_instr_manip_check_trial']



    # Variables that are unique for each stage
    # -------------------------------------------------------------------------

    if stage_id == 'practice':
        pass
    elif stage_id == 'ip_procedure':
        n_staircase_steps = config[stage_id]['n_staircase_trial']

        t_ls = config[stage_id]['t_l']

    elif stage_id == 'experiment':

        def mask(df, key, value):
            return df[df[key] == value]
        pd.DataFrame.mask = mask

        n_reps = config[stage_id]['n_reps']
        m_s_step = config['ip_procedure']['step_size']

        # Read in indifference point procedure log file for extracting
        # indifference points
        ipp_data = pd.read_csv(config['log']['performance']['trial']['file'])

        # Only select data from the last IP procedure iteration; earlier
        # iterations, if any, did not meet preset performance requirements (
        # hence the IP procedure was repeated)
        if ipp_data.empty:
            print('No indifference points available for subject {0:d} (group: {1:d})'.\
                format(config['subject']['subject_ix'],
                       config['subject']['group_ix']))
            pp.core.quit()

        ipp_data = \
            ipp_data[(ipp_data.block_id.str.startswith('i')) & \
                     (ipp_data.iter_ix == max(ipp_data.iter_ix.unique())) & \
                     (ipp_data.trial_type == 'standard')]

        # Determine indifference points - the experimental stimuli are
        # computed based on them
        indifference_points = \
            ipp_data.\
                groupby(by='t_l')['m_s'].\
                agg('last').\
                to_dict()

        # Also, selected the range of delays used for creating experimental
        # stimuli
        t_ls = sorted(indifference_points.keys())


    # Make core trial list
    # =========================================================================

    # Variables that vary between trials
    # -------------------------------------------------------------------------

    if stage_id == 'practice':

        trial_list = \
            pd.DataFrame(
                list(it.product(['standard'],
                                config[stage_id]['m_s'])),
                columns=['trial_type','m_s'])

        trial_list['m_s_cat'] = 'NA'
        trial_list['m_l'] = config[stage_id]['m_l']
        trial_list['t_s'] = config[stage_id]['t_s']
        trial_list['t_l'] = config[stage_id]['t_l']

    elif stage_id == 'ip_procedure':
        n_trial_per_t_l = \
            n_staircase_steps + \
            n_catch_ss + \
            n_catch_ll + \
            n_instr_check

        trial_type = n_catch_ss * ['catch_ss'] + \
                     n_catch_ll * ['catch_ll'] + \
                     n_instr_check * ['instr_check'] + \
                     n_staircase_steps * ['standard']

        # Delay and trial_type vary between trials
        trial_list = \
            pd.DataFrame(
                list(it.product([None]*n_trial_per_t_l,
                                t_ls)),
                columns=['trial_type', 't_l'])

        trial_list['m_s_cat'] = 'NA'

        # Shuffle trial order to make catch trials unpredictable
        for t_l_ix, t_l in enumerate(t_ls):
            random.shuffle(trial_type)
            trial_list.loc[trial_list['t_l'] == t_l, 'trial_type'] = trial_type

    elif stage_id == 'experiment':

        # Shuffle order of frames (blocks), so it varies randomly between
        # participants
        random.shuffle(time_frames)

        # We will present three types of SS amounts: below, at, and above IP
        m_s_cats = ['below_ip', 'at_ip', 'above_ip']

        trial_list = \
            pd.DataFrame(
                list(it.product(time_frames,
                                m_s_cats,
                                range(n_reps),
                                t_ls)),
                columns=['framing', 'm_s_cat', 'i_rep', 't_l'])

        trial_list['m_s'] = None

        # Add m_s
        for i_t_l, t_l in enumerate(t_ls):

            adjustment_factor = 2 * (i_t_l + 1)
            m_s = {'below_ip': indifference_points[t_l] - adjustment_factor * m_s_step,
                   'at_ip': indifference_points[t_l],
                   'above_ip': indifference_points[t_l] + adjustment_factor * m_s_step
                   }

            # If a participant in the indifference procedure always chooses
            # the SS option or always chooses the LL option, then the values
            #  of 'below_ip' and 'above_ip', respectively, should be
            # adjusted, because they will be equal to 0.00 (always chooses
            # SS option) or to the LL amount (always chooses LL option). In
            # these cases, we set m_s one step-size away from the
            # indifference point:

            if m_s['below_ip'] == 0:
                adjustment_factor = (i_t_l + 1)
                m_s['below_ip'] = indifference_points[t_l] - adjustment_factor * m_s_step
            elif m_s['above_ip'] == config['ip_procedure']['m_l']:
                adjustment_factor = (i_t_l + 1)
                m_s['above_ip'] = indifference_points[t_l] + adjustment_factor * m_s_step

            for m_s_cat in m_s_cats:

                trial_list.loc[(trial_list.t_l == t_l) & \
                               (trial_list.m_s_cat == m_s_cat),
                                'm_s'] = m_s[m_s_cat]

    # Variables that remain constant between trials
    # -------------------------------------------------------------------------
    trial_list['session_ix'] = config['session']['session_ix']
    trial_list['waited_for_trigger'] = False
    trial_list['trial_ix'] = None # Will be filled in below
    trial_list['m_unit'] = config['ip_procedure']['m_unit']
    trial_list['t_unit'] = config['ip_procedure']['t_unit']

    if stage_id in ['practice', 'ip_procedure']:
        time_frames = 'neutral'
        trial_list['framing'] = time_frames
        block_ix = 0
        trial_list['block_ix'] = block_ix
        trial_list['block_id'] = '%s%.3d' % (stage_id, block_ix)
    elif isinstance(time_frames, str):
        trial_list['framing'] = time_frames
        block_ix = 0
        trial_list['block_ix'] = block_ix
        trial_list['block_id'] = '%s%.3d' % (stage_id, block_ix)
    else:
        # Add block_ix, block_id
        for block_ix, frame in enumerate(time_frames):
            trial_list. \
                loc[(trial_list.framing == frame), 'block_ix'] = \
                block_ix
            trial_list. \
                loc[(trial_list.framing == frame), 'block_id'] = \
                '%s%.3d' % (stage_id, block_ix)


    if stage_id == 'practice':
        trial_list['t_s'] = config[stage_id]['t_s']
        trial_list['framing'] = 'neutral'

    elif stage_id == 'ip_procedure':

        block_ix = 0
        trial_list['block_ix'] = block_ix
        trial_list['block_id'] = '%s%.3d' % (stage_id, block_ix)
        trial_list['m_s'] = None
        trial_list['m_l'] = config[stage_id]['m_l']


        trial_list['framing'] = 'neutral'

        trial_list['t_s'] = config[stage_id]['t_s']


    elif stage_id == 'experiment':
        trial_list['m_l'] = config['ip_procedure']['m_l']
        trial_list['t_s'] = config['ip_procedure']['t_s']

    # Add catch trials and instruction manipulation check trials, if any
    # ---------------------------------------------------------------------
    if stage_id == 'experiment':
        trial_list['trial_type'] = 'standard'
        if n_catch_ss > 0:
            catch_ss_trials = trial_list.sample(n=n_catch_ss)
            catch_ss_trials.trial_type = 'catch_ss'
            # options have identical amounts, but differ in delay:
            # participants should logically prefer the SS option
            catch_ss_trials.m_s = catch_ss_trials.m_l
            catch_ss_trials.i_rep = None
        else:
            catch_ss_trials = trial_list.sample(n=0)

        if n_catch_ll > 0:
            catch_ll_trials = trial_list.sample(n=n_catch_ll)
            catch_ll_trials.trial_type = 'catch_ll'
            # SS option amount equals zero:
            # participants should logically prefer the LL option
            catch_ll_trials.m_s = 0
            catch_ll_trials.i_rep = None
        else:
            catch_ll_trials = trial_list.sample(n=0)

        if n_instr_check > 0:
            instr_check_trials = trial_list.sample(n=n_instr_check)
            instr_check_trials.trial_type = 'instr_check'
            instr_check_trials.i_rep = None
        else:
            instr_check_trials = trial_list.sample(n=0)

        trial_list = trial_list.append(catch_ss_trials)
        trial_list = trial_list.append(catch_ll_trials)
        trial_list = trial_list.append(instr_check_trials)

    # Make sure that certain variables are represented as ints rather than
    # floats
    trial_list['block_ix'] = trial_list['block_ix'].astype(int)
    trial_list['t_s'] = trial_list['t_s'].astype(int)
    trial_list['t_l'] = trial_list['t_l'].astype(int)



    # Sort rows and columns
    # =========================================================================

    # Sort rows
    # -------------------------------------------------------------------------


    if stage_id in ['practice', 'ip_procedure']:
        trial_list = \
            trial_list. \
                sort(columns=['t_l'], axis='index'). \
                reset_index(drop=True)
    elif stage_id == 'experiment':
        # Randomize all trials, then sort by block index. This ensure that
        # within blocks, trial order is random
        trial_list = \
            trial_list. \
                sample(frac=1). \
                sort(columns=['block_ix'], axis='index'). \
                reset_index(drop=True)

    # Fill in trial_ix
    trial_list['trial_ix'] = range(trial_list.shape[0])
    trial_list['trial_ix'] = trial_list['trial_ix'].astype(int)

    # Sort columns
    # -------------------------------------------------------------------------

    if stage_id in ['practice', 'ip_procedure']:
        col_order = ['session_ix', 'block_ix', 'block_id', 'trial_ix',
                     'waited_for_trigger', 'framing',
                     'trial_type', 'm_s', 'm_s_cat', 'm_l', 'm_unit', 't_s',
                     't_l', 't_unit']
    elif stage_id == 'experiment':
        col_order = ['session_ix', 'block_ix', 'block_id', 'trial_ix',
                     'waited_for_trigger', 'framing',
                     'trial_type', 'm_s', 'm_s_cat', 'm_l', 'm_unit', 't_s',
                     't_l', 't_unit']

    trial_list = trial_list[col_order]

    return trial_list
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
        if block_type == 'practice' or block_type == 'experiment' or \
                block_type == 'break':
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
    trigger_keys = kb['settings']['trigger_keys']
    esc_keys = kb['settings']['esc_keys']
    instruction_stim = config['stimuli']['instruction'][block_type]

    stim_ix = 0


    if block_type == 'block_repeat':
        stim_list = [instruction_stim_ix]
    elif block_type == 'start' or \
            block_type == 'practice' or \
            block_type == 'experiment' or \
            block_type == 'end' or \
            block_type == 'ip_procedure' or \
            block_type == 'break':
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

        if block_type == 'end':
            pp.core.wait(15)
            return

        while no_key_pressed:

            # Collect responses
            rd_key_count, other_keys_pressed = \
                collect_response(rd, kb, other_keys=toggle_keys + trigger_keys)

            # If user pressed key move to next stimulus
            # if sum(rd_key_count.values()) > 0:
            #     window.flip(clearBuffer=True)
            #     stim_ix += 1
            #     break

            if other_keys_pressed:

                if stim_ix < len(stim_list) -1:

                    if other_keys_pressed in toggle_keys:
                        window.flip(clearBuffer=True)
                        if other_keys_pressed == toggle_keys[0]:
                            stim_ix -= 1
                            if stim_ix < 0:
                                stim_ix = 0
                            break
                        elif other_keys_pressed == toggle_keys[1]:
                            stim_ix += 1
                            break
                    else:
                        break
                else:
                    if other_keys_pressed in [toggle_keys[0]] + trigger_keys:
                        window.flip(clearBuffer=True)
                        if other_keys_pressed == toggle_keys[0]:
                            stim_ix -= 1
                            if stim_ix < 0:
                                stim_ix = 0
                            break
                        elif other_keys_pressed in trigger_keys:
                            return
                    else:
                        break



def run_block(config,block_id,trial_list,block_log):
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

    Returns
    -------
    block_log            : pandas.core.frame.DataFrame
                        block-level performance log

    all_crit_met          : bool
                        whether or not all performance criteria have been met

    """

    # 1. PROCESS INPUTS
    # 1.1. Get task configurations
    # 1.2. Define block and trial variables
    # 1.3. Define indifference point procedure variables
    # 2. RUN TRIALS
    # 2.1. Define trial-specific variables
    # 2.1.#. m_s
    # 2.1.#. wait_for_trigger
    # 2.1.#. side of choice options
    # 2.1.#. others
    # 2.2. Present stimuli and collect responses
    # 2.3. Log timings
    # 3. COMPUTE BLOCK STATISTICS
    # 4. RETURN VARIABLES

    # <editor-fold desc="1. PROCESS INPUTS">
    # 1. PROCESS INPUTS
    ###########################################################################
    # <editor-fold desc="1.1. Get task configurations">
    # 1.1. Get task configurations
    # =========================================================================
    frame_rate = config['window']['frame_rate']
    stimuli = config['stimuli']
    ss_opt_pos = config['stim_config']['ss_opt']['pos']
    ll_opt_pos = config['stim_config']['ll_opt']['pos']
    ss_rect_pos = config['stim_config']['ss_rect']['pos']
    ll_rect_pos = config['stim_config']['ll_rect']['pos']


    trial_stats = config['statistics']['trial']
    ip_procedure = config['ip_procedure']
    trial_eval_data = config['evaluation']['trial']
    min_feedback_dur = config['feedback']['trial']['min_duration']
    sess_columns = config['log']['performance']['sess_columns']
    sess_data = config['log']['performance']['sess_data']
    trial_cols = config['log']['performance']['trial']['columns']

    if block_id.startswith('p'):
        stage_id = 'practice'
    elif block_id.startswith('i'):
        stage_id = 'ip_procedure'
    elif block_id.startswith('e'):
        stage_id = 'experiment'

    self_paced = config[stage_id]['self_paced']
    isi_dur = config[stage_id]['isi_dur']
    iti_dur = config[stage_id]['iti_dur']
    min_rt = config[stage_id]['min_rt']
    max_rt = config[stage_id]['max_rt']

    # </editor-fold> # desc="1.1. Get task configurations"
    # <editor-fold desc="1.2. Define block and trial variables">
    # 1.2. Define block and trial variables
    # =========================================================================

    # Block index
    block_ix = trial_list.iloc[0]['block_ix']

    # Block onset: -Inf means that the block starts immediately
    block_ons = -float('inf')

    # Mock times, so that trial starts immediately
    trial_timing = {'ons': -float('inf'),
                    'dur': 0,
                    'self_paced': self_paced,
                    'isi_dur': isi_dur,
                    'min_rt': min_rt,
                    'max_rt': max_rt,
                    'min_feedback_dur': config['feedback']['trial'][
                        'min_duration'],
                    'iti_dur': iti_dur,
                    'refresh_time': 1 / frame_rate}
    t_trial_end = -np.inf

    # Put trial_list indices in a list for easy reference
    trial_list_ixs = trial_list.index.tolist()

    # Init trial_log
    trial_log = pd.DataFrame(index=trial_list_ixs,
                             columns=trial_cols)

    # </editor-fold desc="1.2. Define block and trial variables"> #
    # <editor-fold desc="1.3. Define indifference point procedure variables">
    # 1.3. Define indifference point procedure variables
    # =========================================================================


    if block_id.startswith('i'):

        i_step = iter(range(ip_procedure['n_staircase_trial']))
        m_l = ip_procedure['m_l']

        step_size = m_l * 2 ** -(next(i_step) + 1)
        adjustment_factor = -1
        m_s_staircase = m_l + adjustment_factor * step_size

    # </editor-fold desc="1.3. Define indifference point procedure variables">
    # </editor-fold desc="1. PROCESS INPUTS"> #
    # <editor-fold desc="2. RUN TRIALS">
    # 2. RUN TRIALS
    ###########################################################################
    for trial_list_ix in trial_list_ixs:

        # <editor-fold desc="2.1. Define trial-specific variables">
        # 2.1. Define trial-specific variables
        # =====================================================================
        # <editor-fold desc="2.1.1. Initiate this_trial_log">
        # 2.1.1. Initiate variable for logging this trial's events
        # ---------------------------------------------------------------------
        this_trial_log = \
            copy_series_values(var_names=['session_ix', 'block_ix', 'iter_ix',
                                          'trial_ix',
                                          'trial_type', 'waited_for_trigger',
                                          'framing', 'm_s_cat', 'm_l',
                                          'm_unit', 't_s',
                                          't_l', 't_unit'],
                               source_series=trial_list.loc[trial_list_ix],
                               target_series=pd.Series(index=trial_cols))
        this_trial_log['block_id'] = block_id
        this_trial_log[sess_columns] = sess_data

        # </editor-fold desc="2.1.1. Initiate this_trial_log">
        # <editor-fold desc="2.1.2. Set side of choice stimulus">
        # 2.1.2. Set side of choice stimulus
        # ---------------------------------------------------------------------

        # Get default positions of smaller-sooner and larger-later options
        x_ss_opt, y_ss_opt = ss_opt_pos
        x_ll_opt, y_ll_opt = ll_opt_pos
        x_ss_rect, y_ss_rect = ss_rect_pos
        x_ll_rect, y_ll_rect = ll_rect_pos

        # Randomly determine side, then update and log position
        if random.random() < 0.5:
            stimuli['ss_opt'][0].pos = (x_ss_opt, y_ss_opt)
            stimuli['ll_opt'][0].pos = (x_ll_opt, y_ll_opt)
            stimuli['ss_rect'][0].pos = (x_ss_rect, y_ss_rect)
            stimuli['ll_rect'][0].pos = (x_ll_rect, y_ll_rect)
            this_trial_log['ll_side'] = 'right'
        else:
            # Swap sides
            stimuli['ss_opt'][0].pos = (x_ll_opt, y_ll_opt)
            stimuli['ll_opt'][0].pos = (x_ss_opt, y_ss_opt)
            stimuli['ss_rect'][0].pos = (x_ll_rect, y_ll_rect)
            stimuli['ll_rect'][0].pos = (x_ss_rect, y_ss_rect)
            this_trial_log['ll_side'] = 'left'

        # </editor-fold desc="2.1.2. Set side of choice stimulus">
        # <editor-fold desc="2.1.3. Determine monetary amount SS option">
        # 2.1.3. Determine monetary amount of the smaller-sooner option
        # ---------------------------------------------------------------------
        if block_id.startswith('i'):

            trial_type = trial_list.loc[trial_list_ix, 'trial_type']

            if trial_type in ['catch_ss', 'catch_ll', 'instr_check']:

                # Set monetary amount of SS option
                m_s = {'catch_ss': ip_procedure['m_l'],
                       'catch_ll': 0,
                       'instr_check': m_s_staircase}[trial_type]

            elif trial_type == 'standard':
                m_s = m_s_staircase

        else:
            m_s = trial_list.loc[trial_list_ix, 'm_s']

        this_trial_log['m_s'] = m_s

        # </editor-fold desc="2.1.3. Determine monetary amount SS option">
        # <editor-fold desc="2.1.4. Set stimulus text">
        # 2.1.4. Set stimulus text
        # ---------------------------------------------------------------------
        for stimulus in ['choice_instr', 'ss_opt', 'll_opt']:
            stimuli[stimulus][0].setText(
                make_stim_str(stim=stimulus,
                              framing=trial_list.loc[
                                  trial_list_ix, 'framing'],
                              m_unit=trial_list.loc[
                                  trial_list_ix, 'm_unit'],
                              m_s=m_s,
                              m_l=trial_list.loc[trial_list_ix, 'm_l'],
                              t_unit=trial_list.loc[
                                  trial_list_ix, 't_unit'],
                              t_s=trial_list.loc[trial_list_ix, 't_s'],
                              t_l=trial_list.loc[trial_list_ix, 't_l'])
            )
        # </editor-fold desc="2.1.4. Set stimulus text">

        # </editor-fold> # desc="2.1. Define trial-specific variables"
        # <editor-fold> # desc="2.2. Present stimuli and collect responses"
        # 2.2. Present stimuli and collect responses
        # =====================================================================
        pp.core.wait(t_trial_end + trial_timing['iti_dur'] - pp.core.getTime())

        this_trial_log = \
            run_trial(config=config,
                      window=config['window']['window'],
                      hub=config['apparatus']['hub'],
                      rd=config['apparatus']['rd'],
                      kb=config['apparatus']['kb'],
                      wait_for_trigger=trial_list.loc[trial_list_ix,
                                                      'waited_for_trigger'],
                      trial_log=this_trial_log,
                      trial_timing=trial_timing,
                      trial_stats=trial_stats,
                      trial_eval_data=trial_eval_data,
                      stimuli=stimuli,
                      min_feedback_dur=min_feedback_dur
                      )

        t_trial_end = pp.core.getTime()

        # </editor-fold> # desc="2.2. Present stimuli and collect responses"
        # <editor-fold> # desc="2.3. Log timings"
        # 2.3. Log timings
        # =====================================================================

        # Session timing
        sm, ss = divmod(this_trial_log['trial_ons'], 60)
        sh, sm = divmod(sm, 60)
        this_trial_log['t_session'] = '%d:%02d:%02d' % (sh, sm, ss)

        # Block timing
        if trial_list_ix == trial_list_ixs[0]:
            block_ons = this_trial_log['trial_ons']
        bs, bms = divmod(this_trial_log['trial_ons'] - block_ons,1)
        bm, bs = divmod(bs, 60)
        bh, bm = divmod(bm, 60)
        this_trial_log['t_block'] = '%d:%02d:%02d.%03d' % (bh, bm, bs,bms*1000)

        # </editor-fold> # desc="2.3. Log timings"
        # <editor-fold> # desc="2.4. Save this_trial_log"
        # 2.4. Save this_trial_log
        # =====================================================================

        # Put trial data into data frame and file
        # ---------------------------------------------------------------------
        trial_log.loc[trial_list_ix] = this_trial_log



        with open(config['log']['performance']['trial']['file'],'a+') as fileObj:
            pd.DataFrame(this_trial_log).T.to_csv(fileObj,
                                                  index=False,
                                                  header=False,
                                                  na_rep=np.nan)


        # </editor-fold> # desc="2.4. Save this_trial_log"
        # <editor-fold> # desc="2.5. Update monetary amount of SS option"
        # 2.5. Update monetary amount of smaller-sooner option
        # =====================================================================

        if block_id.startswith('i'):

            # Only update m_s and adjustment factor when trial was a
            # standard trial
            if trial_type in ['catch_ss', 'catch_ll', 'instr_check']:
                pass
            elif trial_type == 'standard':

                try:
                    step_size = m_l * 2 ** -(next(i_step) + 1)

                    adjustment_factor = \
                        {'ss': -1,
                         'll': 1}.get(this_trial_log['choice'], 0)
                except StopIteration:
                    i_step = iter(range(ip_procedure['n_staircase_trial']))
                    m_l = ip_procedure['m_l']

                    step_size = m_l * 2 ** -(next(i_step) + 1)
                    adjustment_factor = -1
                    m_s_staircase = m_l

                # Set adjustment factor according ot choice; if no choice is
                # made, do not adjust m_s

                # If LL chosen on this trial, then the SS offer will be
                # increased on the next trial; if SS chosen on this trial,
                # then the SS offer will be decreased on the next trial. If
                # no choice is made, do not adjust m_s and repeat this trial.

                m_s_staircase = m_s_staircase + adjustment_factor * step_size

                if __debug__:
                    print('Step size: {0:s} {1:.2f}'.
                          format(ip_procedure['m_unit'], step_size))
                    print('Adjustment factor: {0:d}'.
                          format(adjustment_factor))
                    print('Monetary amount SS option next trial: {0:s}{1:.2f}'.
                          format(ip_procedure['m_unit'], + m_s_staircase))
        else:
            pass


        # </editor-fold> # desc="2.5. Update monetary amount of SS option"
    # </editor-fold desc="2. RUN TRIALS">
    # <editor-fold desc="3. COMPUTE BLOCK STATISTICS">
    # 3. COMPUTE BLOCK STATISTICS

    # Compute block stats
    # =========================================================================
    if block_id.startswith('i'):
        df = trial_log[trial_log.block_id == block_id]
        all_crit_met = evaluate_block(config,
                                      df=df,
                                      block_id = block_id,
                                      block_log = block_log)
    else:
        all_crit_met = True

    ###########################################################################
    # </editor-fold desc="3. COMPUTE BLOCK STATISTICS">
    # <editor-fold desc="4. RETURN VARIABLES">
    # 4. RETURN VARIABLES
    return block_log, all_crit_met
    # </editor-fold desc="4. RETURN VARIABLES">
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
                'ip_procedure' or 'experiment'

    trial_list   : pandas.core.frame.DataFrame
                list of trials to present

    """

    block_ixs = trial_list['block_ix'].unique().astype(int)
    block_cols = config['log']['performance']['block']['columns']
    block_log = pd.DataFrame(index=block_ixs,
                             columns=block_cols)

    sess_columns = config['log']['performance']['sess_columns']
    sess_data = config['log']['performance']['sess_data']
    performance_req = config['performance_requirements'][stage_id]

    for block_ix in block_ixs:

        all_crit_met = False
        iter_ix = 0
        force_repeat = performance_req['force_repeat']
        max_n_iter = performance_req['max_n_iter']

        block_id = '%s%.3d' % (stage_id[0], block_ix)

        this_block_log = pd.DataFrame(index=[block_ix], columns=block_cols)

        this_block_log.loc[block_ix, sess_columns] = sess_data
        this_block_log.loc[block_ix, 'block_id'] = block_id
        this_block_log.loc[block_ix, 'block_ix'] = block_ix
        this_block_log.loc[block_ix, 'iter_ix'] = iter_ix

        while not all_crit_met:

            print 'Run %s block %d' % (stage_id, block_ix)

            trial_list_block = trial_list[trial_list['block_ix'] == block_ix]
            trial_list_block.loc[trial_list_block['block_ix'] == block_ix,
                                 'iter_ix'] = iter_ix

            trial_list_block['iter_ix'] = trial_list_block[
                                              'iter_ix'].astype(int)

            # present_instruction(config, stage_id, block_ix)

            this_block_log, all_crit_met = \
                run_block(config=config,
                          block_id=block_id,
                          trial_list=trial_list_block,
                          block_log=this_block_log)

            # Write block log
            # ---------------------------------------------------------
            with open(config['log']['performance']['block']['file'],
                      'a+') as fileObj:
                this_block_log.to_csv(fileObj, index=False, header=False,
                                      na_rep=np.nan)

            if (stage_id == 'experiment') & (block_ix < block_ixs[-1]):
                present_instruction(config, 'break', block_ix)

            if force_repeat:
                if not all_crit_met:
                    if iter_ix == (max_n_iter - 1):
                        present_instruction(config, 'block_repeat', 1)
                        pp.core.wait(15)
                        pp.core.quit()
                    else:
                        iter_ix = iter_ix + 1

                        # Warn subject that block will be repeated
                        present_instruction(config, 'block_repeat', 0)

                        if config['instruction']['enable']:
                            # Present instruction
                            present_instruction(config, stage_id, 0)

            else:
                break
def run_trial(config,wait_for_trigger, hub,trial_log,
              trial_timing,window,rd,kb,trial_stats,trial_eval_data,min_feedback_dur,stimuli):
    """
    Runs a trial

    Parameters
    ----------
    config          : dict
                    specifies StPy experiment properties

    wait_for_trigger  : bool
                    whether or not trial onsent is contingent on trigger

    hub             : psychopy.iohub.client.ioHubConnection
                    interface to the ioHub

    trial_log        : pandas.core.frame.Series
                    empty trial log

    trial_timing     : dict
                    specifies timing properties, such as trial onset, trial
                    duration, inter-trial interval duration, and refresh times

    window          : psychopy.visual.window.Window
                    PsychoPy window object, in which stimuli are presented

    rd              : dict
                    specifies response device properties

    kb              : dict
                    specifies keyboard properties

    trial_stats      : dict
                    specifies which descriptive statistics need to be computed

    trial_eval_data   : dict
                    specifies all data necessary to evaluate trials

    min_feedback_dur     : int or float
                    minimum duration of trial feedback (in seconds)

    stimuli         : dict
                    all stimuli used in the experiment

    Returns
    -------
    trial_log        : pandas.core.frame.Series
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

    # TODO: Rename this to trial_onset?
    t_start = pp.core.getTime()

    # Clear events
    # -------------------------------------------------------------------------
    hub.clearEvents('all')

    if __debug__:
        t_events_cleared = pp.core.getTime()
        print '* Events cleared: t = %f ms; dt = %f ms' % \
              (1000 * (t_events_cleared - t_start),
               1000 * (t_events_cleared - t_start))

    # Present stimuli
    # -------------------------------------------------------------------------

    if trial_log['trial_type'] == 'instr_check':
        default_instruction = stimuli['choice_instr'][0].text
        manip_instruction = default_instruction.replace('Choose between:',
                                                      'Press return:')
        stimuli['choice_instr'][0].setText(manip_instruction)


    # Draw choice instruction
    stimuli['choice_instr'][0].draw()

    # Initially, present choice instruction only to let participants process it
    choice_instr_ons = window.flip()

    # Draw choice options, too
    stimuli['choice_instr'][0].draw()
    stimuli['ss_opt'][0].draw()
    stimuli['ll_opt'][0].draw()

    pp.core.wait(trial_timing['isi_dur'])

    # Next, show choice options too
    opt_ons = window.flip()

    # Log
    trial_log['trial_ons'] = choice_instr_ons
    trial_log['trial_dur'] = 'NA'

    # Collect responses
    # -------------------------------------------------------------------------

    trial_log = collect_response(rd=rd, kb=kb, log=trial_log,
                                 t0=opt_ons,
                                 min_rt=trial_timing['min_rt'],
                                 max_rt=trial_timing['max_rt'])

    # Set feedback duration
    if trial_log['rt'] < trial_timing['min_rt']:
        feedback_dur = trial_timing['min_feedback_dur'] + \
                       trial_timing['min_rt'] - \
                       trial_log['rt']
    else:
        feedback_dur = trial_timing['min_feedback_dur']

    # Evaluate trial
    # -------------------------------------------------------------------------
    trial_log = evaluate_trial(eval_data=trial_eval_data,
                               feedback_dur=feedback_dur,
                               window=window,
                               stimuli=stimuli,
                               log=trial_log)

    # Present intertrial interval
    # -------------------------------------------------------------------------
    # stimuli['iti'][0].setAutoDraw(True)
    stimuli['iti'][0].setText(' ')
    stimuli['iti'][0].draw()
    window.flip()

    # Wrap up
    # -------------------------------------------------------------------------
    return trial_log