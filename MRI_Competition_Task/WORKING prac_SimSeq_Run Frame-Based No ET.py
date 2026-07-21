from psychopy import prefs, plugins, sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware, monitors
from psychopy.constants import (NOT_STARTED, STARTED, STOPPED, FINISHED, PRESSED, RELEASED, FOREVER, priority)
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy.hardware import keyboard
import math
from math import sin, cos, pi, radians
import numpy as np
from collections import deque
from string import ascii_letters, digits
import csv
import random
import itertools
import os
import json
import sys
import time
import pylink
from psychopy.tools.monitorunittools import pix2deg
logging.console.setLevel(logging.ERROR)

###### PARAMETERS ######################################################################################################################

# Initialize the global clock and keyboard
globalClock = core.Clock()
kb = keyboard.Keyboard(clock = globalClock)

# Experiment
ATTENTION_CONDS = ['FIX', 'COV']
BLOCK_DESIGN = [('RVF','SIM'),('LVF','SEQ')]
NUM_RUNS = 6 # how many runs in one feature condition
RUNS_PER_COND = int(NUM_RUNS//len(ATTENTION_CONDS)) # equal number of FIX and COV runs per feature condition (runs 1,2,3 are FIX; 4,5,6 are COV)
NUM_SIM_BLOCKS = 1 # per run
NUM_SEQ_BLOCKS = 1 # per run
NUM_BLANK_BLOCKS = 2 # one before and one after each run
NUM_TRIALS = 3 # per block

# Timing
BLANK_BLOCK_DURATION = 16 # seconds
PERIPH_STIM_DURATION = 1 # seconds
RSVP_RATE = 0.5 #sec; durations of RSVP pokemon presentation
TRIAL_DURATION = PERIPH_STIM_DURATION*4 # sec
PSTIM_TARGET_FREQ = [1,3] # pstim color targets will occur every 1-3 trials (4-12s)
POKEMON_TARGET_FREQ = [7,15] # pokemon targets will occur every 15-30 pokemon (3.75-7.5s) in the RSVP
RESPONSE_WINDOW = 1.5 # sec; responses during this window this will be coded as hits

# Stim parameters
DISTANCE = 60 # cm, pt distance from screen
PERIPHERAL_STIM_SIZE = 1.25 #DVA; size of each peripheral stimulus (circle)
POKEMON_SIZE = [1.5, 1.5] # DVA, size of the pokemon during RSVP
POKEMON_POS = (0,0) # location of rsvp pokemon
NUM_PSTIMS = 4 # number of peripheral stimuli
GRID_SIZE = 4 #DVA; height and width of the peripheral stimulus grid
ECCENTRICITY = 7 #DVA from the center of the screen to the center of the grid
PERIPHERAL_STIM_COLORS = {"red":(0.8027, 0.4268, 0.6013), "orange":(0.6965,0.5192,0.2497), "green":(0.4837,0.6019,0.3269), 
    "cyan":(0.2239,0.6331,0.6725), "blue":(0.4307, 0.5730, 0.9193), "magenta":(0.7312, 0.4429, 0.8911)} # from CIELUV space
CLR_SPC = 'rgb1'
FREQUENCY =  1.0 # Hz
AMPLITUDE = 0.75 # half of the stim's total motion in DVA
ANGLES = [0, 30, 60, 90, 120, 150] # all possible angles of motion
TARGET_ANGLE = 90 # vertical motion
TARGET_COLOR = 'red' # (0.8027, 0.4268, 0.6013)
GAZE_BOUND = 3 # if gaze shifts more than this from fixation point, flag the trial

# Response key
RESPONSE_KEY = '1'
SCANNER_KEY = 'equal' 

###### SETUP ###########################################################################################################################

# Options to choose from in the dialogue box
pokemon_names = ["Bulbasaur", "Pikachu", "Squirtle", "Charmander", "Magikarp", "Raticate", "Pidgey",
    "Metapod", "Jigglypuff", "Butterfree", "Psyduck", "Caterpie", "Krabby",
    "Haunter", "Vulpix", "Eevee", "Sandshrew", "Wartortle", "Charmeleon", "Clefairy",
    "Ponyta", "Mankey"]
    
feat_conds = ['color', 'motion', 'color-motion']

# Present dialogue box to collect experiment parameters
exp_name = 'SimSeq'
exp_info = {
    'Participant ID': '9999', 
    'Session': '001',
    'Pokemon': pokemon_names,
    'Condition': feat_conds,
    'Runs': 'FIX: 1,2,3; COV: 4,5,6'
}

dlg = gui.DlgFromDict(dictionary=exp_info, title=exp_name)
if dlg.OK == False:
    core.quit()
    sys.exit()

# Get and print inputted parameters
target_pokemon = exp_info['Pokemon'].strip().capitalize() 
feat_cond = exp_info['Condition']
run_list = [int(x.strip()) for x in exp_info['Runs'].split(',')] # which runs to display
sub_id = exp_info['Participant ID']

print('###############################################################################\n')
print(f'Performing runs {run_list} of {feat_cond} condition.')
print('Target pokemon:', target_pokemon)
print('Target color:', TARGET_COLOR)
print('\n###############################################################################\n')

# Establish data output directory
time_str = time.strftime("%m_%d_%Y", time.localtime())
root_dir = os.path.dirname(os.path.abspath(__file__))
participant_folder = os.path.join(root_dir, 'data', f"{exp_info['Participant ID']}_{exp_name}_Session{exp_info['Session']}_{time_str}")
os.makedirs(participant_folder, exist_ok=True)
feature_folder = os.path.join(participant_folder, f"{feat_cond}_cond_{target_pokemon}")
os.makedirs(feature_folder, exist_ok=True)
trials_filename = os.path.join(feature_folder, f"exp_trials_{exp_info['Participant ID']}_Session{exp_info['Session']}_{feat_cond}_{target_pokemon}")
rsvp_filename = os.path.join(feature_folder, f"exp_rsvp_{exp_info['Participant ID']}_Session{exp_info['Session']}_{feat_cond}_{target_pokemon}")

# Create experiment handlers to manage the data files
thisExp = data.ExperimentHandler(name=exp_name, version='', extraInfo=exp_info,
                                runtimeInfo=None, originPath=os.path.abspath(__file__),
                                savePickle=True, saveWideText=True,
                                dataFileName=trials_filename)

rsvpExp = data.ExperimentHandler(name='rsvp',extraInfo=exp_info,
                                savePickle=True, saveWideText=True,
                                dataFileName=rsvp_filename)

# Set column orders to improve data file readability
trials_column_order = ['welcome.start','instructions.start', 'instructions.end', 'ready_screen.start', 'ready_screen.end', 'block.start', 'block.end', 'feat_cond', 'run', 'attention_cond','block',
    'presentation_cond', 'vf', 'trial', 'rsvp_seq', 'pstim_colors', 'pstim_angles', 'pstim_phases','trial.start', 'trial.end', 'pstim.onset', 
    'target_shown', 'target.onset', 'press_times', 'rts', 'keypresses', 'hit']
    
if feat_cond == 'color':
    pstim_order = ['blue_None.onset', 'blue_None.offset', 'cyan_None.onset', 'cyan_None.offset', 'green_None.onset', 'green_None.offset', 
                'magenta_None.onset', 'magenta_None.offset', 'orange_None.onset', 
                'orange_None.offset', 'red_None.onset', 'red_None.offset']
elif feat_cond == 'motion':
    pstim_order = ['black_0.onset', 'black_0.offset', 'black_30.onset', 'black_30.offset', 'black_60.onset', 'black_60.offset', 'black_90.onset', 
                'black_90.offset', 'black_120.onset', 'black_120.offset', 'black_150.onset', 'black_150.offset']
elif feat_cond == 'color-motion':
    pstim_order = ['blue_0.onset', 'blue_0.offset', 'cyan_0.onset', 'cyan_0.offset', 'magenta_0.onset', 'magenta_0.offset', 'orange_0.onset', 'orange_0.offset', 'green_0.onset', 'green_0.offset', 
             'red_0.onset', 'red_0.offset', 'red_30.onset', 'red_30.offset',
             'blue_30.onset', 'blue_30.offset', 'cyan_30.onset', 'cyan_30.offset', 'green_30.onset' , 'green_30.offset','magenta_30.onset', 'magenta_30.offset', 'orange_30.onset', 'orange_30.offset', 
             'cyan_60.onset', 'cyan_60.offset', 'green_60.onset', 'green_60.offset', 'blue_60.onset' , 'blue_60.offset', 'magenta_60.onset', 'magenta_60.offset', 'orange_60.onset', 'orange_60.offset', 
             'red_60.onset', 'red_60.offset','red_90.onset', 'red_90.offset', 'blue_120.onset', 'blue_120.offset', 'cyan_120.onset', 'cyan_120.offset', 'green_120.onset', 'green_120.offset', 
             'magenta_120.onset', 'magenta_120.offset', 'orange_120.onset', 'orange_120.offset', 'blue_150.onset', 'blue_150.offset', 'cyan_150.onset', 
             'cyan_150.offset', 'green_150.onset', 'green_150.offset', 'magenta_150.onset', 'magenta_150.offset', 'orange_150.onset', 'orange_150.offset', 
             'red_150.onset', 'red_150.offset']

trials_cols = trials_column_order + pstim_order
for col in trials_cols:
    thisExp.addData(col, '')
    
rsvp_column_order = ['ready_screen.start', 'ready_screen.end', 'block.start', 'rsvp_seq', 'run', 'block', 'trial']

for col in rsvp_column_order:
    rsvpExp.addData(col, '')
    
# Mark the experiment as started
exp_info['expDate'] = data.getDateStr(format = '%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6)
thisExp.status = STARTED
rsvpExp.status = STARTED

# Window setup (will need to be adjusted to match the display monitor)
Eizo = monitors.Monitor('Eizo', width = 51.84, distance = DISTANCE) # screen width (cm) and distance from the screen (cm)
Eizo.setSizePix([1920, 1200])
win = visual.Window(fullscr=True, color=[0.9032,0.8051,0.9655],
            size=Eizo.getSizePix(), screen=1,
            winType='pyglet', allowStencil=False,
            monitor=Eizo, colorSpace=CLR_SPC,
            backgroundImage='', backgroundFit='none',
            blendMode='avg', useFBO=False,
            units='deg', 
            checkTiming=False)

# Get monitor's refresh rate (CRITICAL FOR PROPER RSVP PRESENTATION)
frame_rate = win.getActualFrameRate()
if frame_rate is None:
    frame_rate = 60.0 # fallback to 60Hz
    print("Could not measure frame rate, defaulting to 60Hz")
exp_info['frameRate'] = frame_rate

# Load Pokémon images from folder
image_dir = os.path.join(root_dir, 'images')
pokemon_dict = {name: visual.ImageStim(win, name=name, image=os.path.join(image_dir, f"{i+1:03}.png"))
    for i, name in enumerate(pokemon_names)} # dictionary of the pokemon where the key is their name and the values are the ImageStim object

# Calculate the coordinates of the center of the peripheral grid based on the display monitor
cent2cent_spacing = GRID_SIZE - PERIPHERAL_STIM_SIZE # distance from center to center of peripheral stimuli
offset = cent2cent_spacing / 2 # how much to move in x and y from the center of the grid to the center of each peripheral stimulus
angle_rad = np.deg2rad(45) # polar angle from x axis to the center of the grid in radians
gridcent_x = ECCENTRICITY * np.cos(angle_rad) # X coordinate of the grid center in RVF, make negative for LVF
gridcent_y = ECCENTRICITY * np.sin(angle_rad) # Y coordinate of the grid center

# Calculate the peripheral stimuli positions in the 2x2 grid relative to the center of the grid calculated above
rvf_topleft = [gridcent_x - offset, gridcent_y + offset]  # Top left peripheral stimulus
rvf_topright = [gridcent_x + offset, gridcent_y + offset]  # Top right peripheral stimulus
rvf_botleft= [gridcent_x - offset, gridcent_y - offset]  # Bottom left peripheral stimulus
rvf_botright = [gridcent_x + offset, gridcent_y - offset]  # Bottom right peripheral stimulus
lvf_topleft = [-gridcent_x - offset, gridcent_y + offset]  # Top left peripheral stimulus in LVF
lvf_topright = [-gridcent_x + offset, gridcent_y + offset]  # Top right peripheral stimulus in LVF
lvf_botleft = [-gridcent_x - offset, gridcent_y - offset]  # Bottom left peripheral stimulus in LVF
lvf_botright = [-gridcent_x + offset, gridcent_y - offset]  # Bottom right peripheral stimulus in LVF


###### FUNCTIONS #######################################################################################################################

def terminate_task():
    """ 
    Helper function to save data and close the window.
    """
    
    # Mark end of experiment
    print("Ending experiment.")
    thisExp.nextEntry()
    thisExp.addData('experiment.end', globalClock.getTime(format='float'))
    rsvpExp.nextEntry()
    rsvpExp.addData('experiment.end', globalClock.getTime(format='float'))
    
    # Save trial and rsvp data
    thisExp.saveAsWideText(trials_filename + '.csv', delim='auto')
    thisExp.saveAsPickle(trials_filename)
    rsvpExp.saveAsWideText(rsvp_filename + '.csv', delim='auto')
    rsvpExp.saveAsPickle(rsvp_filename)
    logging.flush()
    print(f"Data saved to {os.path.abspath(feature_folder)}")
    
    # Clear psychopy window
    if win is not None:
        win.clearAutoDraw()
        win.flip()
        
    # Change experiment status to finished
    thisExp.status = FINISHED
    rsvpExp.status = FINISHED
    thisExp.abort()
    rsvpExp.abort()

    # Close the window
    win.close()
    core.quit()
    sys.exit()
    
def generate_rsvps():
    """" 
    Generates lists of pokemon names to be presented during each run. 
    
    Returns a dictionary.
    Keys are run numbers. Value for each key is a list of pokemon names that serve as the RSVP stimulus 
    presentation list for that entire run. 
    
    Checks if dictionaries have already been generated for this subID, given feature condition
    and target pokemon are the same. 
    
    Prevents back to back presentations of the same pokemon. 
    """
    
    # Check if dictionaries exist in feature folder. If so, load and return them
    filename = os.path.join(feature_folder, f"{sub_id}_rsvps.json")
    
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        rsvps = {int(k): v for k, v in data['exp_rsvps'].items()}
        print(f"Loaded existing RSVPs from {filename}")
        return rsvps
    
    # Initialize dictionaries and calculate number of pokemon needed
    rsvps = {}
    run_dur = (NUM_BLANK_BLOCKS*BLANK_BLOCK_DURATION) + (((NUM_SIM_BLOCKS+NUM_SEQ_BLOCKS)*NUM_TRIALS)*TRIAL_DURATION)
    pokemon_per_run = int(run_dur // RSVP_RATE) # total number of pokemon in a run
    
    # Create all of the lists of pokemon names for experiment runs. 
    for run_idx in range(NUM_RUNS//2):
        # Get times for target occurrences in this block
        target_indices = []
        next_target_idx = random.randint(*POKEMON_TARGET_FREQ)
        while next_target_idx < pokemon_per_run:
            target_indices.append(next_target_idx)
            next_target_idx += random.randint(*POKEMON_TARGET_FREQ)
        # Build sequence, prevent back-to-back repeats
        pokemon_list = []
        for idx in range(pokemon_per_run):
            if idx in target_indices:
                pokemon_list.append(target_pokemon)
            else:
                distractor_options = [p for p in pokemon_names if p != target_pokemon]
                if idx > 0:
                    distractor_options = [p for p in distractor_options if p != pokemon_list[-1]] #ensure no back to back pokemon
                pokemon_list.append(random.choice(distractor_options))
        rsvps[run_idx+1] = pokemon_list
  
    # RSVPs for runs 4,5,6 are the same as for 1,2,3 - EDIT TO MAKE IT MORE ROBUST
    rsvps[4] = rsvps[1]
    rsvps[5] = rsvps[2]
    rsvps[6] = rsvps[3]
    
    # Save RSVPs to file for future loading
    data = {
        'exp_rsvps': {str(k): v for k, v in rsvps.items()}
    }
    with open(filename, 'w') as f:
        json.dump(data, f)
    print(f"Generated and saved RSVPs to {filename}")
    
    return rsvps
    
def generate_sim_onsets():
    """
    Generates pseudo-random onset times for peripheral stimulus grid onset in SIM trials. 
    
    Returns a dictionary. 
    Keys are run numbers. Value for each key is a list of the onset times for all SIM trials of that run.
    
    Example: sim_onsets[1] is list of integers, one for each SIM trial of the first run. 
    
    Checks if dictionaries have already been generated for this subID, given feature condition
    and target pokemon are the same. 
    
    Onsets will never be back to back (i.e., if the grid was presented during the last second of the previous trial,
    it will not be presented during the first second of the current trial). 
    """
    
    # Check if dictionaries exist in feature folder. If so, load and return them
    filename = os.path.join(feature_folder, f"{sub_id}_sim_onsets.json")
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        sim_onsets = {int(k): v for k, v in data['exp_sim_onsets'].items()}
        print(f"Loaded existing sim onsets from {filename}")
        return sim_onsets
    
    # Generate a list of integers (0-3) with length = total number of SIM trials. The integers
    # are seconds within the trial when the peripheral grid can be presented.
    sim_onsets = {}
    for run in range(1, RUNS_PER_COND+1):
        trial_onsets = []
        for trial_idx in range(NUM_SIM_BLOCKS*NUM_TRIALS):
            if trial_idx != 0 and trial_onsets[trial_idx-1] == 3:
                trial_onsets.append(random.choice([1, 2, 3]))
            else:
                trial_onsets.append(random.choice([0, 1, 2, 3]))
        sim_onsets[run] = trial_onsets
        
    # SIM onsets for runs 1,2,3 are the same for 4,5,6 - EDIT TO MAKE IT MORE ROBUST
    sim_onsets[4] = sim_onsets[1]
    sim_onsets[5] = sim_onsets[2]
    sim_onsets[6] = sim_onsets[3]
    
    # Save SIM onsets to file for future reloading
    data = {
        'exp_sim_onsets': {str(k): v for k, v in sim_onsets.items()}
    }
    with open(filename, 'w') as f:
        json.dump(data, f)
    print(f"Generated and saved sim onsets to {filename}")
    
    return sim_onsets
    
def assign_grids(feat_cond):
    """ 
    Returns a dictionary. 
    Keys are run numbers. Value for each key are lists of 
        
    Example:
        Calling trial_grids[1][32] will give the grid layout (list of colors and angles) for 
        the 33rd trial in the first run. 
        
    Checks if dictionaries have already been generated for this subID, given feature condition
    and target pokemon are the same. 
    """
    
    # Check if dictionaries exist in feature folder. If so, load and return them
    filename = os.path.join(feature_folder, f"{sub_id}_grids.json")
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        trial_grids = {int(k): v for k, v in data['exp_trial_grids'].items()}
        print(f"Loaded existing grids from {filename}\n")
        print('###############################################################################\n')
        return trial_grids

    # Initialize dictionaries and calculate how many trials are in the run
    trial_grids = {}
    total_trials = NUM_TRIALS * (NUM_SEQ_BLOCKS + NUM_SIM_BLOCKS)

    #-------Color feature condition
    all_color_combos = list(itertools.combinations(PERIPHERAL_STIM_COLORS.keys(), NUM_PSTIMS)) # every possible combination of 4 colors from list of all possible colors
    all_color_configs = []
    for combo in all_color_combos:
        all_color_configs.extend(itertools.permutations(combo)) # every possible ordering of those 4 colors
    random.shuffle(all_color_configs)
    
    # Separate the list of all possible permutations (configs)
    # NOTE: target color will only appear in target trials and in the same position in the grid (closest to fixation)
    rvf_target_colors = [grid for grid in all_color_configs if TARGET_COLOR in grid and grid[2] == TARGET_COLOR] # target color in RVF target position
    lvf_target_colors = [grid for grid in all_color_configs if TARGET_COLOR in grid and grid[3] == TARGET_COLOR] # target color in LVF target position
    colors_without_target = [grid for grid in all_color_configs if TARGET_COLOR not in grid] # no target color
    
    #-------Motion and color-motion condition
    all_angle_combos = list(itertools.combinations(ANGLES, NUM_PSTIMS)) # every possible combination of 4 angles from list of all possible angles
    all_angle_configs = []
    for combo in all_angle_combos:
        all_angle_configs.extend(itertools.permutations(combo)) # every possible ordering of those 4 angles
    random.shuffle(all_angle_configs)
    
    # Separate the list of all possible permutations (configs)
    # NOTE: target angle will only appear in target trials and in the same position in the grid (closest to fixation)
    rvf_target_angles = [c for c in all_angle_configs if c[2] == TARGET_ANGLE] # target angle in RVF target position
    lvf_target_angles = [c for c in all_angle_configs if c[3] == TARGET_ANGLE] # target angle in LVF target position
    angles_without_target = [c for c in all_angle_configs if TARGET_ANGLE not in c] # no target angle

    #-------Helper function to assign grids for a single run
    def assign_run_grids():
        nonlocal total_trials
        
        # Generate indices for trials that should have the target based on the target frequency
        target_trials = []
        target_trial_idx = random.randint(*PSTIM_TARGET_FREQ)
        while target_trial_idx <= total_trials:
            target_trials.append(target_trial_idx)
            target_trial_idx += random.randint(*PSTIM_TARGET_FREQ)

        run_assignments = []
        used_color_configs = set()
        used_angle_configs = set()

        # Assign a grid to each trial based on its VF and whether the trial idx is in target_trials
        for trial_idx in range(total_trials):
            block_num = (trial_idx // NUM_TRIALS) % len(BLOCK_DESIGN)
            block_vf = BLOCK_DESIGN[block_num][0]
            is_target = trial_idx in target_trials
            chosen_colors = None
            chosen_angles = None

            if feat_cond == 'color':
                if is_target:
                    available = [c for c in (rvf_target_colors if block_vf == 'RVF' else lvf_target_colors) if c not in used_color_configs] or (rvf_target_colors if block_vf == 'RVF' else lvf_target_colors)
                else:
                    available = [c for c in colors_without_target if c not in used_color_configs] or colors_without_target
                chosen_colors = random.choice(available)
                used_color_configs.add(chosen_colors)
                run_assignments.append({'colors': chosen_colors, 'angles': [None] * NUM_PSTIMS})

            elif feat_cond == 'motion':
                if is_target:
                    available = [a for a in (rvf_target_angles if block_vf == 'RVF' else lvf_target_angles) if a not in used_angle_configs] or (rvf_target_angles if block_vf == 'RVF' else lvf_target_angles)
                else:
                    available = [a for a in angles_without_target if a not in used_angle_configs] or angles_without_target
                chosen_angles = random.choice(available)
                used_angle_configs.add(chosen_angles)
                run_assignments.append({'colors': ('black', 'black', 'black', 'black'), 'angles': chosen_angles})

            elif feat_cond == 'color-motion':
                if is_target:
                    color_available = [c for c in (rvf_target_colors if block_vf == 'RVF' else lvf_target_colors) if c not in used_color_configs] or (rvf_target_colors if block_vf == 'RVF' else lvf_target_colors)
                    angle_available = [a for a in (rvf_target_angles if block_vf == 'RVF' else lvf_target_angles) if a not in used_angle_configs] or (rvf_target_angles if block_vf == 'RVF' else lvf_target_angles)
                else:
                    # NOTE: Non-target trials may show the target color in the target position (closest to fixation)
                    # but will never have the target angle
                    color_pool = colors_without_target + (rvf_target_colors if block_vf == 'RVF' else lvf_target_colors)
                    color_available = [c for c in color_pool if c not in used_color_configs] or color_pool
                    angle_available = [a for a in angles_without_target if a not in used_angle_configs] or angles_without_target
                chosen_colors = random.choice(color_available)
                chosen_angles = random.choice(angle_available)
                used_color_configs.add(chosen_colors)
                used_angle_configs.add(chosen_angles)
                run_assignments.append({'colors': chosen_colors, 'angles': chosen_angles})

        return run_assignments

    #-------Assign grids for experiment runs
    for run in range(1, RUNS_PER_COND + 1):
        trial_grids[run] = assign_run_grids()

    # Runs 4,5,6 reuse the same grids as runs 1,2,3 - EDIT TO MAKE IT MORE ROBUST
    trial_grids[4] = trial_grids[1]
    trial_grids[5] = trial_grids[2]
    trial_grids[6] = trial_grids[3]

    # Save to file for future reloading
    data = {
        'exp_trial_grids': {str(k): v for k, v in trial_grids.items()}
    }
    with open(filename, 'w') as f:
        json.dump(data, f)
    print(f"Generated and saved grids to {filename}\n")
    print('###############################################################################\n')

    return trial_grids

def create_trial_dicts(rsvp, run_sim_onsets, trial_grids):
    """
    Call at the start of each run. 
    
    Returns a list of trial dictionaries for all of the trials in a run. Helps capture the
    trial variables in the trials data file.
    
    Args are the generated rsvp, sim onsets, and trial grids for the entire run. 
        
    Returns:
        trial_dicts (list): Each value is a dictionary of trial vraibles, including:
            'block_num', 'trial_num', 'visual_field', 'presentation_cond', 
            'peripheral_grid', 'rsvp_seq'
            
    Example:
        Calling trial_dicts[32] within a run, returns the trial dictionary of
        the 33rd trial in that run.
    """    
    
    trial_dicts = []

    # Divide the rsvp into sub-lists of pokemon to be presented each trial
    pokemon_per_trial = int(TRIAL_DURATION//RSVP_RATE)
    pokemon_per_blank = int(BLANK_BLOCK_DURATION//RSVP_RATE)
    trial_pokemon = rsvp[pokemon_per_blank:-pokemon_per_blank] # remove pokemon presented during blank blocks
    trial_rsvps = [trial_pokemon[i:i+pokemon_per_trial] for i in range(0, len(trial_pokemon), pokemon_per_trial)] # divide rsvp into sublists, one per trial
    
    # Convert sim_onsets from a list to a deque to use .popleft
    sim_onsets = deque(run_sim_onsets)
        
    for block_idx, (visual_field, present_cond) in enumerate(BLOCK_DESIGN):
        for trial_in_block in range(NUM_TRIALS):
            trial_idx = block_idx * NUM_TRIALS + trial_in_block
            trial_dict = {
                'block_num': block_idx + 1,
                'trial_num': trial_idx + 1,
                'visual_field': visual_field,
                'presentation_cond': present_cond,
                'peripheral_grid': trial_grids[trial_idx],
                'rsvp_seq': trial_rsvps[trial_idx]
            }
            if present_cond == 'SIM':
                trial_dict['grid_onset'] = sim_onsets.popleft()
            trial_dicts.append(trial_dict)
    return trial_dicts

def show_instructions(feat_cond, attention_cond):
    """
    Helper function to display instructions based on feature condition and attention condition.
    Waits until space or escape is pressed. 
    """
    
    thisExp.addData('instructions.start', globalClock.getTime(format='float'))
    
    # COV instructions vary depending on the feature condition
    if feat_cond == 'color':
        cov_instructions_text = visual.TextStim(win=win, text=(
                "There's a Pokémon Party happening right now, and the Pokémon are hungry!\n\n\n\n\n\n"
                f"The Pokémon like to eat {TARGET_COLOR} circles like these! Can you help feed them?\n\n\n"
                f"Press the button as fast as you can every time you see a {TARGET_COLOR} circle.\n\n\n"
                "Ready to start playing?"), font='Arial', units='deg', pos=(0, 0), height=1, wrapWidth=1400, 
                color='black', colorSpace=CLR_SPC)
        instruc_pstim = visual.Circle(win, pos = rvf_botleft, radius = PERIPHERAL_STIM_SIZE/2, 
            units = 'deg', anchor = 'center', fillColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], lineColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], colorSpace = CLR_SPC)
        instruc_pstim2 = visual.Circle(win, pos = lvf_botright, radius = PERIPHERAL_STIM_SIZE/2, 
            units = 'deg', anchor = 'center', fillColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], lineColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], colorSpace = CLR_SPC)
    elif feat_cond == 'motion':   
        cov_instructions_text = visual.TextStim(win=win, text=(
                "There's a Pokémon Party happening right now, and the Pokémon are hungry!\n\n\n\n\n\n"
                f"The Pokémon like to eat black circles that move up and down like this! Can you help feed them?\n\n\n"
                f"Press the button as fast as you can every time you see a black circle moving up and down.\n\n\n"
                "Ready to start playing?"), font='Arial', units='deg', pos=(0, 0), height=1, wrapWidth=1400, 
                color='black', colorSpace=CLR_SPC)
        instruc_pstim = visual.Circle(win, pos = rvf_botleft, radius = PERIPHERAL_STIM_SIZE/2, 
            units = 'deg', anchor = 'center', fillColor='black', lineColor='black')
        instruc_pstim2 = visual.Circle(win, pos = lvf_botright, radius = PERIPHERAL_STIM_SIZE/2, 
            units = 'deg', anchor = 'center', fillColor='black', lineColor='black')
    elif feat_cond == 'color-motion':
        cov_instructions_text = visual.TextStim(win=win, text=(
                "There's a Pokémon Party happening right now, and the Pokémon are hungry!\n\n\n\n\n\n"
                f"The Pokémon like to eat {TARGET_COLOR} circles that move up and down like this! Can you help feed them?\n\n\n"
                f"Press the button as fast as you can every time you see a {TARGET_COLOR} circle moving up and down.\n\n\n"
                "Ready to start playing?"), font='Arial', units='deg', pos=(0, 0), height=1, wrapWidth=1400, 
                color='black', colorSpace=CLR_SPC)
        instruc_pstim = visual.Circle(win, pos = rvf_botleft, radius = PERIPHERAL_STIM_SIZE/2, 
            units = 'deg', anchor = 'center', fillColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], lineColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], colorSpace = CLR_SPC)
        instruc_pstim2 = visual.Circle(win, pos = lvf_botright, radius = PERIPHERAL_STIM_SIZE/2, 
            units = 'deg', anchor = 'center', fillColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], lineColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], colorSpace = CLR_SPC)

    # FIX instructions are the same across all feature conditions
    fix_instructions_text = visual.TextStim(win=win, text=(
        "There's a Pokémon Party happening right now, and the Pokémon are playing hide and seek!\n\n"
        f"The Pokémon are trying to find {target_pokemon}! Can you help them? {target_pokemon} will show up like this:\n\n\n"
        f"\n\nPress the button as fast as you can every time you see {target_pokemon}.\n\n\n"
        "Ready to start playing?"), font='Arial', units='deg', pos=(0, 0), height=1, wrapWidth=1400, 
        color='black', colorSpace=CLR_SPC)
    
    # Draw instruction components on the screen, wait until space or escape is pressed
    if attention_cond == "FIX":
        fix_instructions_text.draw()
        pokemon_dict[target_pokemon].pos = (0, 0) 
        pokemon_dict[target_pokemon].size = [1.5, 1.5]
        pokemon_dict[target_pokemon].draw()
        win.flip()
        keys = event.waitKeys(keyList=['space', 'escape'])
        
        if 'escape' in keys:
            terminate_task()
        
    elif attention_cond == "COV":
        if feat_cond != 'color': # dots will be moving if feat_cond is not color
            base_pos = instruc_pstim.pos
            base_pos2 = instruc_pstim2.pos
            while True:
                t = globalClock.getTime()
                dy = AMPLITUDE * np.sin(2 * np.pi * FREQUENCY * t)
                instruc_pstim.pos = base_pos + np.array([0, dy])
                instruc_pstim2.pos = base_pos2 + np.array([0, dy])
                cov_instructions_text.draw()
                instruc_pstim.draw()
                instruc_pstim2.draw()
                win.flip()
                
                keys = event.getKeys(keyList=['space', 'escape'])
                
                if 'space' in keys:
                    break
                elif 'escape' in keys:
                    terminate_task()
        else:
            cov_instructions_text.draw()
            instruc_pstim.draw()
            instruc_pstim2.draw()
            win.flip()
            keys = event.waitKeys(keyList=['space', 'escape'])

    thisExp.addData('instructions.end', globalClock.getTime(format='float')) 
    
def end_run_screen(feat_cond, attention_cond):
    """ 
    Helper function to display text at the end of a run based on feature condition and attention condition.
    Waits until space or escape is pressed.
    """
    
    thisExp.nextEntry()
    thisExp.addData('end_run_screen.start', globalClock.getTime(format='float'))
    rsvpExp.addData('end_run_screen.start', globalClock.getTime(format='float'))
    
    end_text = visual.TextStim(win=win, text=(), font='Arial', units='deg', pos=(0, 0), height=1.2, wrapWidth=1700, 
        color='black', colorSpace=CLR_SPC)
        
    # FIX end screen is the same across feature conditions
    if attention_cond == 'FIX':
        end_text.text = (f"Great job finding {target_pokemon}!")
        pokemon_dict[target_pokemon].pos = (0, -5) 
        pokemon_dict[target_pokemon].size = (5,5)
        end_text.draw()
        pokemon_dict[target_pokemon].draw()
        win.flip()
        
        keys = event.waitKeys(keyList=['space', 'escape'])
        
        if 'escape' in keys:
            terminate_task()
        
        # Reset size and position of the target pokemon
        pokemon_dict[target_pokemon].size = POKEMON_SIZE
        pokemon_dict[target_pokemon].pos = (0,0)
        
    # COV end screen varies depending on feature condition
    elif attention_cond == 'COV':
        if feat_cond == 'color':
            target_pstim = visual.Circle(win, name = TARGET_COLOR, pos = (0,-5), radius = 1 , units = 'deg', anchor = 'center', 
                fillColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], lineColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], colorSpace = CLR_SPC)
        elif feat_cond == 'motion':
            target_pstim = visual.Circle(win, name = TARGET_COLOR, pos = (0,-5), radius = 1 , units = 'deg', 
                anchor = 'center', fillColor='black', lineColor='black')
        elif feat_cond == 'color-motion':
            target_pstim = visual.Circle(win, name = TARGET_COLOR, pos = (0,-5), radius = 1 , units = 'deg', anchor = 'center', 
                fillColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], lineColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], colorSpace = CLR_SPC)
        end_text.text = (f"Great job feeding the Pokémon the circles!")
        
        if feat_cond == 'motion' or feat_cond == 'color-motion':
            base_pos = target_pstim.pos
            while True:
                t = globalClock.getTime()
                dy = AMPLITUDE * np.sin(2 * np.pi * FREQUENCY * t)
                target_pstim.pos = base_pos + np.array([0, dy])
                end_text.draw()
                target_pstim.draw()
                win.flip()
                
                keys = event.getKeys(keyList=['space', 'escape'])
                
                if 'space' in keys:
                    break
                elif 'escape' in keys:
                    terminate_task()
        else:
            end_text.draw()
            target_pstim.draw()
            win.flip()
            
            keys = event.waitKeys(keyList=['space', 'escape'])
            
            if 'escape' in keys:
                terminate_task()
                
    thisExp.addData('end_run_screen.end', globalClock.getTime(format='float'))
    rsvpExp.addData('end_run_screen.end', globalClock.getTime(format='float'))
    
def ready_screen():
    """ Text screen before the start of the experiment trials."""
    
    thisExp.addData('ready_screen.start', globalClock.getTime(format='float'))
    rsvpExp.addData('ready_screen.start', globalClock.getTime(format='float'))

    ready_text = visual.TextStim(win=win, text=("Ready to start the real game?"), font='Arial', units='deg', pos=(0, 0), height=1.2, wrapWidth=1700, 
        color='black', colorSpace=CLR_SPC)
        
    ready_text.draw()
    win.flip()
    
    keys = event.waitKeys(keyList=[SCANNER_KEY,'escape'])
    if 'escape' in keys:
        terminate_task()
        
    thisExp.addData('ready_screen.end', globalClock.getTime(format='float'))
    rsvpExp.addData('ready_screen.end', globalClock.getTime(format='float'))

#def run_trial(feat_cond, attention_cond, run, rsvp, rsvp_state, trial_dict, last_target_onset):    
#    """
#    Function to run a single trial.
#    """
#    
#    # Reset variables for trial
#    kb.clearEvents()
#    pstim_idx = 0
#    hit = 0
#    press_times = []
#    rts = []
#    target_onset_recorded = False
#    flip_time = None
#    
#    # Get presentation condition
#    is_seq_trial = trial_dict['presentation_cond'] == 'SEQ'
#    is_sim_trial = trial_dict['presentation_cond'] == 'SIM'
#    
#    # Set up rsvp sequence variables
#    frames_per_pokemon = round(RSVP_RATE*frame_rate)
#    pokemon_per_trial = round(TRIAL_DURATION/RSVP_RATE)
#    total_trial_frames = pokemon_per_trial * frames_per_pokemon
#    rsvp_sequence = trial_dict['rsvp_seq']
#    
#    # Extract peripheral stim locations
#    vf = trial_dict['visual_field']
#    if vf[0] == 'R':
#        grid_positions = [rvf_topright, rvf_topleft, rvf_botleft, rvf_botright]
#    elif vf[0] == 'L':
#        grid_positions = [lvf_topright, lvf_topleft, lvf_botleft, lvf_botright]
#
#    # Map peripheral stim colors, angles, and phases to their positions 
#    pstim_to_draw = []
#    pstim_info = {}
#    pstim_colors = trial_dict['peripheral_grid']['colors']
#    pstim_angles = trial_dict['peripheral_grid']['angles']
#    phases = [0, 0, np.pi/2, np.pi/2] # half of the circles will start at middle of motion, and half will start at a peak
#    random.shuffle(phases)
#    
#    for color, angle, pos, phase in zip(pstim_colors, pstim_angles, grid_positions, phases):
#        if color == 'black':
#            pstim = visual.Circle(win, name = f"{color}_{angle}", pos = pos, radius = PERIPHERAL_STIM_SIZE/2 , 
#                units = 'deg', anchor = 'center', fillColor=color, lineColor=color, colorSpace=CLR_SPC)
#        else:
#            pstim = visual.Circle(win, name = f"{color}_{angle}", pos = pos, radius = PERIPHERAL_STIM_SIZE/2 , 
#                units = 'deg', anchor = 'center', fillColor=PERIPHERAL_STIM_COLORS[color], lineColor=PERIPHERAL_STIM_COLORS[color], colorSpace=CLR_SPC)
#        pstim_to_draw.append(pstim) # create a list of pstims to draw in this trial
#        pstim_info[pstim.name] = {'angle': angle, 'base_pos': pstim.pos, 'phase': phase}
#    pstim_onset_recorded_dict = {pstim.name: False for pstim in pstim_to_draw}
#    pstim_offset_recorded_dict = {pstim.name: False for pstim in pstim_to_draw}
#        
#    # Save trial variables to the data file
#    thisExp.addData('feat_cond', feat_cond)
#    thisExp.addData('run', f'exp{run}')
#    thisExp.addData('block', trial_dict['block_num'] )
#    thisExp.addData('trial', trial_dict['trial_num'])
#    thisExp.addData('attention_cond', attention_cond)
#    thisExp.addData('presentation_cond', trial_dict['presentation_cond'])
#    thisExp.addData('vf', vf[0])
#    thisExp.addData('rsvp_seq', rsvp_sequence)
#    thisExp.addData('pstim_colors', pstim_colors)
#    thisExp.addData('pstim_angles', pstim_angles)
#    thisExp.addData('pstim_phases', phases)
#
#    # Set peripheral stim grid onsets
#    if is_sim_trial:
#        pstim_start_time = trial_dict['grid_onset'] # start time for entire grid in sim trials
#        thisExp.addData('pstim.onset', pstim_start_time)
#    elif is_seq_trial:
#        pstim_onsets = np.arange(0, TRIAL_DURATION, PERIPH_STIM_DURATION).tolist() # seconds; SEQ condition: each peripheral stimuli will be shown one at a time for PERIPH STIM DURATION (s)
#        thisExp.addData('pstim.onset', pstim_onsets[0]) # get the first pstim onset
#    
#    # Check whether target is present this trial
#    if attention_cond == 'FIX': # target shown if target pokemon appears in rsvp during the trial
#        target_shown = target_pokemon in rsvp_sequence 
#    elif attention_cond == 'COV':
#        if feat_cond == 'color': # target shown if circle closest to fix is target color
#            if vf[0] == 'L':
#                target_shown = TARGET_COLOR == pstim_colors[3]
#            elif vf[0] == 'R':
#                target_shown = TARGET_COLOR == pstim_colors[2]
#        elif feat_cond == 'motion': # target shown if black circle closest to fix is moving vertically
#            if vf[0] == 'L':
#                target_shown = TARGET_ANGLE == pstim_angles[3]
#            elif vf[0] == 'R':
#                target_shown = TARGET_ANGLE == pstim_angles[2]
#        elif feat_cond == 'color-motion': # target shown if circle closest to fix is moving vertically and is target color
#            if vf[0] == 'L':
#                target_shown = (TARGET_ANGLE == pstim_angles[3] and TARGET_COLOR == pstim_colors[3])
#            elif vf[0] == 'R':
#                target_shown = (TARGET_ANGLE == pstim_angles[2] and TARGET_COLOR == pstim_colors[2])
#                
#    # Blank block if this is the first trial in the run
#    if trial_dict['trial_num'] == 1:
#        last_target_onset = run_blank_block(feat_cond, attention_cond, run, 1, rsvp, rsvp_state, last_target_onset)
#
#    while rsvp_state['frame'] % frames_per_pokemon != 0:
#        if rsvp_state['current_name'] is not None:
#            current_pokemon = pokemon_dict[rsvp_state['current_name']]
#            current_pokemon.draw()
#        win.flip()
#        rsvp_state['frame'] += 1
#
#    trial_start = globalClock.getTime()
#    thisExp.addData('trial.start', trial_start)
#
#    current_pokemon = None
#    
#    # start of trial loop
#    for frameN in range(total_trial_frames):
#        t = globalClock.getTime()
#        
#        # -------------------- RSVP ------------------------------------------------------------------------------
#        record_pokemon_onset = False
#        update_target_onset = False
#        
#        # Update current pokemon at every frames_per_pokemon interval
#        if rsvp_state['frame'] % frames_per_pokemon == 0:
#            current_pokemon = pokemon_dict[rsvp[rsvp_state['idx']]]
#            record_pokemon_onset = True # flag this window flip to record pokemon onset
#            
#            # FIX only: Determine whether pokemon is the target pokemon, if so, flag the next window flip to record target onset
#            if attention_cond == "FIX" and target_shown:
#                if current_pokemon.name == target_pokemon:
#                    update_target_onset = True
#        
#        # Continuously draw the current pokemon at every frame
#        if current_pokemon is not None:
#            current_pokemon.draw()
#            
#        # -------------------- Peripheral Stim Grids ------------------------------------------------------------------------------
#        record_pstim_onset = False
#        record_pstim_offset = False
#        
#        # SIM presentation condition
#        if is_sim_trial:
#            # Continuously update the position of all circles and redraw to make them move
#            for pstim in pstim_to_draw:
#                if trial_start + pstim_start_time <= t < trial_start + pstim_start_time + PERIPH_STIM_DURATION:
#                    angle = pstim_info[pstim.name]['angle']
#                    base_pos = pstim_info[pstim.name]['base_pos']
#                    phase = pstim_info[pstim.name]['phase']
#                    if angle is not None:
#                        elapsed = t - (trial_start + pstim_start_time)
#                        angle_rad = np.deg2rad(angle)
#                        offset = AMPLITUDE * np.sin(2*np.pi*FREQUENCY*elapsed + phase)
#                        dy = offset * np.sin(angle_rad)
#                        dx = offset * np.cos(angle_rad)
#                        pstim.pos = base_pos + (dx,dy)
#                    pstim.draw()
#                    
#                    record_pstim_onset = True
#                    
#                    # COV only: Determine whether target pstim is in the grid, if so, flag the next window flip to record target onset
#                    if attention_cond == "COV" and target_shown:
#                        update_target_onset = True
#                        
#                elif t >= trial_start + pstim_start_time + PERIPH_STIM_DURATION:
#                    record_pstim_offset = True
#            
#        # SEQ presentation condition
#        if is_seq_trial and pstim_idx < len(pstim_to_draw):
#            current_pstim = pstim_to_draw[pstim_idx]
#            onset_time = trial_start + pstim_onsets[pstim_idx]
#            
#            if onset_time <= t < onset_time + PERIPH_STIM_DURATION:
#                # Continuously update the position of circles and redraw to make them move
#                angle = pstim_info[current_pstim.name]['angle']
#                base_pos = pstim_info[current_pstim.name]['base_pos']
#                phase = pstim_info[current_pstim.name]['phase']
#                if angle is not None:
#                    elapsed = t - onset_time
#                    angle_rad = np.deg2rad(angle)
#                    offset = AMPLITUDE * np.sin(2*np.pi*FREQUENCY*elapsed + phase)
#                    dx = offset * np.cos(angle_rad)
#                    dy = offset * np.sin(angle_rad)
#                    current_pstim.pos = base_pos + (dx,dy)
#                current_pstim.draw()
#                
#                record_pstim_onset = True
#                
#                # COV only: Determine whether pstim being drawn is the target, if so, flag the next window flip to record target onset
#                if attention_cond == "COV" and target_shown:
#                    if vf[0] == 'L' and pstim_idx == 3:
#                        update_target_onset = True
#                    elif vf[0] == 'R' and pstim_idx == 2:
#                        update_target_onset = True
#                
#            # After PERIPH_STIM_DURATION (s), move on to next pstim circle
#            elif t >= onset_time + PERIPH_STIM_DURATION:
#                pstim_idx +=1
#                record_pstim_offset = True
#                
#        # --------------------------------Record stimulus onsets and offsets when the window flips ------------------------------------------------------------------------------
#        def on_flip():
#            
#            nonlocal last_target_onset, target_onset_recorded, flip_time
#            
#            flip_time = globalClock.getTime(format='float') # get time when window flipped
#            
#            # Record pstim onsets and offsets for SIM trials, only recorded at first window flip after being drawn/erased
#            if is_sim_trial: 
#                if record_pstim_onset:
#                    for pstim in pstim_to_draw:
#                        if not pstim_onset_recorded_dict[pstim.name]:
#                            thisExp.addData(f'{pstim.name}.onset', flip_time)
#                            pstim_onset_recorded_dict[pstim.name] = True
#                elif record_pstim_offset:
#                    for pstim in pstim_to_draw:
#                        if not pstim_offset_recorded_dict[pstim.name]:
#                            thisExp.addData(f'{pstim.name}.offset', flip_time)
#                            pstim_offset_recorded_dict[pstim.name] = True
#            
#            # Record pstim onsets and offsets for SEQ trials, only recorded at first window flip after being drawn/erased
#            if is_seq_trial and pstim_idx > 0:
#                prev_idx = min(pstim_idx - 1, len(pstim_to_draw) - 1)
#                prev_pstim = pstim_to_draw[prev_idx]
#                
#                if record_pstim_onset and not pstim_onset_recorded_dict[current_pstim.name]:
#                    thisExp.addData(f'{current_pstim.name}.onset', flip_time)
#                    pstim_onset_recorded_dict[current_pstim.name] = True
#                elif record_pstim_offset and not pstim_offset_recorded_dict[prev_pstim.name]:
#                    thisExp.addData(f'{prev_pstim.name}.offset', flip_time)
#                    pstim_offset_recorded_dict[prev_pstim.name] = True
#                    
#            # Record pokemon onset if this is the first window flip that the pokemon is being drawn
#            if record_pokemon_onset:
#                if rsvp_state['current_name'] is not None: # start recording offsets after the first pokemon has been recorded
#                    rsvpExp.addData('stim.offset', flip_time) 
#                    rsvpExp.nextEntry()
#                    
#                # Update rsvp state dictionary with the current pokemon's variables
#                rsvp_state['current_name'] = current_pokemon.name
#                rsvp_state['current_onset'] = flip_time
#                rsvp_state['current_run'] = f'exp{run}'
#                rsvp_state['current_block'] = trial_dict['presentation_cond']
#                rsvp_state['current_trial'] = trial_dict['trial_num']
#                rsvp_state['current_feat_cond'] = feat_cond
#                rsvp_state['current_attention_cond'] = attention_cond
#                
#                rsvpExp.addData('run', rsvp_state['current_run'])
#                rsvpExp.addData('block', rsvp_state['current_block'])
#                rsvpExp.addData('trial', rsvp_state['current_trial'])
#                rsvpExp.addData('feat_cond', rsvp_state['current_feat_cond'])
#                rsvpExp.addData('attention_cond', rsvp_state['current_attention_cond'])
#                rsvpExp.addData('stim', rsvp_state['current_name'])
#                rsvpExp.addData('stim.onset', rsvp_state['current_onset'])
#                
#                rsvp_state['idx'] += 1
#                
#            # Update the last_target_onset if target was drawn during this flip
#            if update_target_onset:
#                last_target_onset = flip_time
#                if attention_cond == 'FIX': # can have multiple targets per trial
#                    rsvpExp.addData('target.onset', last_target_onset)
#                elif attention_cond == 'COV' and not target_onset_recorded: # can only have one target per trial
#                    target_onset_recorded = True
#                thisExp.addData('target.onset', last_target_onset)
#        
#        # Call on_flip function when the window flips and save the last time window was flipped
#        win.callOnFlip(on_flip)
#        win.flip()
#        
#        # --------------------------- Response Check ------------------------------------------------------------------------------
#        # If participant presses response button before RESPONSE_WINDOW seconds have passed since last target onset, count it as a hit
#        keys = kb.getKeys(keyList=[RESPONSE_KEY, 'escape'], waitRelease=False, clear=True)
#        for key in keys:
#            if key.name == 'escape':
#                terminate_task()
#            elif key.name == RESPONSE_KEY:
#                press_time = key.rt
#                press_times.append(press_time)
#                if last_target_onset is not None:
#                    rt = press_time - last_target_onset
#                    rts.append(rt)
#                    if 0 < rt <= RESPONSE_WINDOW:
#                        hit = 1
#                    if attention_cond == "FIX":
#                        rsvpExp.addData('press_time', press_time)
#                        rsvpExp.addData('rt', rt)
#                        rsvpExp.addData('hit', hit)
#                        
#        # Update the rsvp state dictionary
#        rsvp_state['frame'] += 1
#    
#    # ------------------ After the trial loop ------------------------------------------------------------------------------
#    # Save offset for pstims displayed in last 1 second of trial
#    if is_seq_trial and 'current_pstim' in locals():
#        thisExp.addData(f'{current_pstim.name}.offset', flip_time)
#    if is_sim_trial and pstim_start_time == 3:
#        for pstim in pstim_to_draw:
#            thisExp.addData(f'{pstim.name}.offset', flip_time)
#        
#    # Save trial outputs to data file
#    thisExp.addData('target_shown', target_shown)
#    thisExp.addData('press_times', press_times)
#    thisExp.addData('rts', rts)
#    thisExp.addData('hit', hit)
#    thisExp.addData('keypresses', len(press_times)) # how many times pt responded in the trial
#    thisExp.addData('trial.end', flip_time)
#    thisExp.nextEntry()
#        
#    # Print trial data
#    print(f"Run: {run}, Block: {trial_dict['block_num']}, Trial: {trial_dict['trial_num']}, Target shown: {target_shown}, Number of keypresses: {len(press_times)}")
#    
#    # Logic if this was the last trial 
#    last_trial = NUM_TRIALS * (NUM_SEQ_BLOCKS + NUM_SIM_BLOCKS)
#    if trial_dict['trial_num'] == last_trial:
#        # close the final pokemon in the stream 
#        if rsvp_state['current_name'] is not None:
#            rsvpExp.addData('stim.offset', flip_time)
#            rsvpExp.nextEntry()
#            
#            rsvp_state['current_name'] = None
#            rsvp_state['current_onset'] = None
#            rsvp_state['current_run'] = None
#            rsvp_state['current_block'] = None
#            rsvp_state['current_trial'] = None
#            rsvp_state['current_feat_cond'] = None
#            rsvp_state['current_attention_cond'] = None
#            
#        # run the blank block at the end of the run
#        last_target_onset = run_blank_block(feat_cond, attention_cond, run, 2, rsvp, rsvp_state, last_target_onset)
#
#    return last_target_onset
    
#def experiment_run(feat_cond, attention_cond, run, rsvp, sim_onsets, trial_grids):
#    """
#    Function to perform one experiment run. 
#    """
#    # Reset last target onset
#    last_target_onset = None
#    
#    # Create dictionary to keep track of the rsvp stream across trials
#    rsvp_state = {
#        'idx':0, 
#        'frame': 0,
#        'current_name': None,
#        'current_onset': None,
#        'current_run': None,
#        'current_block': None,
#        'current_trial': None,
#        'current_feat_cond': None,
#        'current_attention_cond': None
#    }
#    
#    # Extract visuals for this run
#    trial_dicts = create_trial_dicts(rsvp, sim_onsets, trial_grids) # create trial dicts for this run
#
#    # Run all trials using trial dictionaries, updating accuracy
#    for trial_dict in trial_dicts:
#        last_target_onset = run_trial(feat_cond, attention_cond, run, rsvp, rsvp_state, trial_dict, last_target_onset)
#    
#    # Display text at the end of the run
#    end_run_screen(feat_cond, attention_cond)

def experiment_run(feat_cond, attention_cond, run, rsvp, sim_onsets, trial_grids):
    """
    Function to perform one experiment run.
    """

    # --- Setup ----------------------------------------------------------------------------------------
    frames_per_blank = round(BLANK_BLOCK_DURATION * frame_rate)
    frames_per_trial = round(TRIAL_DURATION * frame_rate)
    frames_per_pokemon = round(RSVP_RATE * frame_rate)
    pokemon_per_blank = int(BLANK_BLOCK_DURATION // RSVP_RATE)

    trial_dicts = create_trial_dicts(rsvp, sim_onsets, trial_grids)
    print(trial_dicts)

    # Map trial onsets/offsets to global run frames (end_frame is exclusive)
    trial_schedule = []
    frame_counter = frames_per_blank     # first trial starts after the blank block
    for trial_dict in trial_dicts:
        trial_schedule.append({
            'trial_dict':  trial_dict,
            'start_frame': frame_counter,
            'end_frame':   frame_counter + frames_per_trial,
        })
        frame_counter += frames_per_trial

    trials_per_run = (NUM_SIM_BLOCKS + NUM_SEQ_BLOCKS) * NUM_TRIALS
    total_run_frames = round((trials_per_run * frames_per_trial) + (NUM_BLANK_BLOCKS * frames_per_blank))
    print('total run frames:', total_run_frames)
    print('end of final blank:', trial_schedule[-1]['end_frame'] + frames_per_blank)

    # --- Trial-level state ----------------------------------------------------------------------------------------
    last_target_onset = None
    target_onset_recorded = False
    trial_offset_recorded = False
    current_trial_idx = 0
    current_trial = trial_schedule[current_trial_idx]
    current_trial_dict = current_trial['trial_dict']
    rsvp_idx = 0
    pstim_clock = core.Clock()

    current_pokemon = None
    is_seq_trial = False
    is_sim_trial = False
    target_shown = False
    rts = []
    press_times = []
    keypresses = 0

    pstim_idx = 0
    pstim_to_draw = []
    pstim_info = {}
    pstim_onset_recorded_dict = {}
    pstim_offset_recorded_dict = {}
    pstim_onsets_frames = []
    pstim_start_frame = None
    current_pstim = None
    vf = None
    pstim_colors = []
    pstim_angles = []
    phases = []

    # --- Helper: shared sinusoidal motion calculation ------------------------------------------------------
    def _apply_motion(pstim, elapsed):
        """Update pstim.pos based on its angle, base_pos, and phase."""
        info = pstim_info[pstim.name]
        angle = info['angle']
        base_pos = info['base_pos']
        phase = info['phase']
        if angle is not None:
            angle_rad = np.deg2rad(angle)
            offset = AMPLITUDE * np.sin(2 * np.pi * FREQUENCY * elapsed + phase)
            pstim.pos = base_pos + (offset * np.cos(angle_rad), offset * np.sin(angle_rad))

    # --- on_flip callback (defined once, outside the frame loop) ──────────────
    # Flags are re-read from the enclosing scope at the moment the flip fires,
    # so defining this outside the loop is safe and avoids recreating it every frame.
    def on_flip():
        nonlocal last_target_onset, rsvp_idx, target_onset_recorded

        flip_time = globalClock.getTime(format='float')

        # Record trial offset/onset
        if record_trial_offset:
            thisExp.addData('trial.end', flip_time)
            thisExp.nextEntry()
        if record_trial_onset:
            pstim_clock.reset() 
            thisExp.addData('trial.start', flip_time)
            thisExp.addData('feat_cond', feat_cond)
            thisExp.addData('run', f'exp{run}')
            thisExp.addData('block', current_trial_dict['block_num'])
            thisExp.addData('trial', current_trial_dict['trial_num'])
            thisExp.addData('attention_cond', attention_cond)
            thisExp.addData('presentation_cond', current_trial_dict['presentation_cond'])
            thisExp.addData('vf', vf[0])
            thisExp.addData('rsvp_seq', current_trial_dict['rsvp_seq'])
            thisExp.addData('pstim.onset', pstim_onset_to_log)
            thisExp.addData('pstim_colors', pstim_colors)
            thisExp.addData('pstim_angles', pstim_angles)
            thisExp.addData('pstim_phases', phases)
            rsvpExp.addData('rsvp_seq', current_trial_dict['rsvp_seq'])
            trial_idx = int(current_trial_dict['trial_num'])-1
            if trial_idx % NUM_TRIALS == 0: # mark the start of each block
                thisExp.addData('block.start', flip_time)
                rsvpExp.addData('block.start', flip_time)

        # Record pokemon onset (and previous pokemon's offset)
        if record_pokemon_onset:
            if rsvp_idx == 0:
                rsvpExp.addData('block.start', flip_time)
            elif rsvp_idx == pokemon_per_blank:
                rsvpExp.addData('block.end', flip_time)
            if rsvp_idx > 0:
                rsvpExp.addData('stim.offset', flip_time)
                rsvpExp.nextEntry()
            rsvpExp.addData('run', f'exp{run}')
            rsvpExp.addData('block', 'blank' if rsvp_idx < pokemon_per_blank or rsvp_idx >= len(rsvp)-pokemon_per_blank else current_trial_dict['block_num'])
            rsvpExp.addData('trial', 0 if rsvp_idx < pokemon_per_blank or rsvp_idx >= len(rsvp)-pokemon_per_blank else current_trial_dict['trial_num'])
            rsvpExp.addData('feat_cond', feat_cond)
            rsvpExp.addData('attention_cond', attention_cond)
            rsvpExp.addData('stim', current_pokemon.name)
            rsvpExp.addData('stim.onset', flip_time)
            rsvp_idx += 1

        # Record pstim onsets/offsets for SIM trials
        # (all pstims share one pstim_start_frame, so onset/offset fires together)
        if is_sim_trial:
            if record_pstim_onset:
                for pstim in pstim_to_draw:
                    if not pstim_onset_recorded_dict[pstim.name]:
                        thisExp.addData(f'{pstim.name}.onset', flip_time)
                        pstim_onset_recorded_dict[pstim.name] = True
            elif record_pstim_offset:
                for pstim in pstim_to_draw:
                    if not pstim_offset_recorded_dict[pstim.name]:
                        thisExp.addData(f'{pstim.name}.offset', flip_time)
                        pstim_offset_recorded_dict[pstim.name] = True

        # Record pstim onsets/offsets for SEQ trials
        if is_seq_trial:
            if pstim_idx > 0:
                prev_pstim = pstim_to_draw[pstim_idx-1]
                if record_pstim_offset and not pstim_offset_recorded_dict[prev_pstim.name]:
                    thisExp.addData(f'{prev_pstim.name}.offset', flip_time)
                    pstim_offset_recorded_dict[prev_pstim.name] = True
            if record_pstim_onset and current_pstim is not None and not pstim_onset_recorded_dict[current_pstim.name]:
                thisExp.addData(f'{current_pstim.name}.onset', flip_time)
                pstim_onset_recorded_dict[current_pstim.name] = True

        # Record target onset
        if update_target_onset:
            last_target_onset = flip_time
            if attention_cond == 'FIX':
                rsvpExp.addData('target.onset', last_target_onset)
            elif attention_cond == 'COV' and not target_onset_recorded:
                target_onset_recorded = True
                thisExp.addData('target.onset', last_target_onset)

    # --- Frame loop ----------------------------------------------------------------------------------------
    for frame in range(total_run_frames):

        # Reset all logging flags
        update_target_onset = False
        record_pokemon_onset = False
        record_pstim_onset = False
        record_pstim_offset = False
        record_trial_onset = False
        record_trial_offset = False

        # --- RSVP ----------------------------------------------------------------------------------------
        if frame % frames_per_pokemon == 0 and rsvp_idx < len(rsvp):
            current_pokemon = pokemon_dict[rsvp[rsvp_idx]]
            record_pokemon_onset = True
            if attention_cond == "FIX" and current_pokemon.name == target_pokemon:
                update_target_onset = True

        if current_pokemon is not None:
            current_pokemon.draw()

        # --- Trial advancement (end_frame is exclusive) ------------------------------------------------------------
        if frame >= current_trial['end_frame']:
            if not trial_offset_recorded:
                record_trial_offset = True
                trial_offset_recorded = True
            if current_trial_idx + 1 < len(trial_schedule):
                current_trial_idx  += 1
                current_trial = trial_schedule[current_trial_idx]
                current_trial_dict = current_trial['trial_dict']

        # --- Trial initialization (runs only on the trial's first frame) --------------------------------------------
        if frame == current_trial['start_frame']:
            record_trial_onset = True
            pstim_idx = 0
            target_onset_recorded = False
            trial_offset_recorded = False
            rts = []
            press_times = []
            keypresses = 0

            is_seq_trial = current_trial_dict['presentation_cond'] == 'SEQ'
            is_sim_trial = current_trial_dict['presentation_cond'] == 'SIM'

            vf = current_trial_dict['visual_field']
            grid_positions = (
                [rvf_topright, rvf_topleft, rvf_botleft, rvf_botright] if vf[0] == 'R'
                else [lvf_topright, lvf_topleft, lvf_botleft, lvf_botright]
            )

            pstim_to_draw = []
            pstim_info = {}
            pstim_colors = current_trial_dict['peripheral_grid']['colors']
            pstim_angles = current_trial_dict['peripheral_grid']['angles']
            phases = [0, 0, np.pi/2, np.pi/2]
            random.shuffle(phases)

            for color, angle, pos, phase in zip(pstim_colors, pstim_angles, grid_positions, phases):
                fill = color if color == 'black' else PERIPHERAL_STIM_COLORS[color]
                pstim = visual.Circle(
                    win, name=f"{color}_{angle}", pos=pos,
                    radius=PERIPHERAL_STIM_SIZE/2, units='deg', anchor='center',
                    fillColor=fill, lineColor=fill, colorSpace=CLR_SPC
                )
                pstim_to_draw.append(pstim)
                pstim_info[pstim.name] = {'angle': angle, 'base_pos': pstim.pos, 'phase': phase}

            pstim_onset_recorded_dict = {p.name: False for p in pstim_to_draw}
            pstim_offset_recorded_dict = {p.name: False for p in pstim_to_draw}

            if is_sim_trial:
                pstim_start_frame = current_trial['start_frame'] + current_trial_dict['grid_onset'] * frame_rate
                pstim_onset_to_log = current_trial_dict['grid_onset']
            elif is_seq_trial:
                pstim_onsets = np.arange(0, TRIAL_DURATION, PERIPH_STIM_DURATION).tolist()
                pstim_onsets_frames = [(onset * frame_rate) + current_trial['start_frame'] for onset in pstim_onsets]
                pstim_onset_to_log = pstim_onsets[0]

            # Check whether target pstim is present this trial
            target_shown = False
            if attention_cond == 'FIX':
                target_shown = target_pokemon in current_trial_dict['rsvp_seq']
            elif attention_cond == 'COV':
                idx = 3 if vf[0] == 'L' else 2   # index of the circle closest to fixation
                if feat_cond == 'color':
                    target_pstim_name = f'{TARGET_COLOR}_None'
                    target_shown = TARGET_COLOR == pstim_colors[idx]
                elif feat_cond == 'motion':
                    target_pstim_name = f'black_{TARGET_ANGLE}'
                    target_shown = TARGET_ANGLE == pstim_angles[idx]
                elif feat_cond == 'color-motion':
                    target_pstim_name = f'{TARGET_COLOR}_{TARGET_ANGLE}'
                    target_shown = (TARGET_COLOR == pstim_colors[idx] and TARGET_ANGLE == pstim_angles[idx])

            print(f"Run: {run}, Block: {current_trial_dict['block_num']}, "
                  f"Trial: {current_trial_dict['trial_num']}, Target shown: {target_shown}")

        # --- SIM: draw all pstims simultaneously --------------------------------------------------------------------------
        if is_sim_trial and pstim_start_frame is not None and frame >= current_trial['start_frame']:
            elapsed = pstim_clock.getTime()
            stim_end = pstim_start_frame + round(PERIPH_STIM_DURATION * frame_rate)
            for pstim in pstim_to_draw:
                if pstim_start_frame <= frame < stim_end:
                    _apply_motion(pstim, elapsed)
                    pstim.draw()
                    record_pstim_onset = True
                    if attention_cond == "COV" and pstim.name == target_pstim_name:
                        if not target_onset_recorded:
                            update_target_onset = True
                elif frame >= stim_end:
                    record_pstim_offset = True  # assumes all pstims share one stim_end

        # --- SEQ: draw one pstim at a time ------------------------------------------------------------------------------
        if is_seq_trial and pstim_idx < len(pstim_to_draw) and frame >= current_trial['start_frame']:
            current_pstim = pstim_to_draw[pstim_idx]
            onset_frame = pstim_onsets_frames[pstim_idx]
            stim_end = onset_frame + round(PERIPH_STIM_DURATION * frame_rate)
            elapsed = pstim_clock.getTime()
            if onset_frame <= frame < stim_end:
                _apply_motion(current_pstim, elapsed)
                current_pstim.draw()
                record_pstim_onset = True
                if attention_cond == "COV" and current_pstim.name == target_pstim_name:
                    if (vf[0] == 'L' and pstim_idx == 3) or (vf[0] == 'R' and pstim_idx == 2):
                        if not target_onset_recorded:
                            update_target_onset = True
            elif frame >= stim_end:
                pstim_idx += 1
                record_pstim_offset = True

        # --- Flip + response check ----------------------------------------------------------------------------------------
        win.callOnFlip(on_flip)
        win.flip()

        keys = kb.getKeys(keyList=[RESPONSE_KEY, 'escape'], waitRelease=False, clear=True)
        for key in keys:
            if key.name == 'escape':
                terminate_task()
            elif key.name == RESPONSE_KEY:
                keypresses += 1
                press_times.append(key.rt)
                if last_target_onset is not None:
                    rt = key.rt - last_target_onset
                    rts.append(rt)
                    print('Key pressed. RT:', rt, rts)
                    if attention_cond == "FIX":
                        rsvpExp.addData('press_time', key.rt)
                        rsvpExp.addData('rt', rt)
                        if 0 < rt <= RESPONSE_WINDOW:
                            rsvpExp.addData('hit', 1)
                    elif attention_cond == "COV":
                        thisExp.addData('rt', rt)
                        thisExp.addData('rts', rts)
                        thisExp.addData('press_times', press_times)
                        thisExp.addData('keypresses', keypresses)
                        if 0 < rt <= RESPONSE_WINDOW:
                            thisExp.addData('hit', 1)

    end_run_screen(feat_cond, attention_cond)
  
###### WELCOME SCREEN ##################################################################################################################

# Clear the window and set start time
win.flip() 
thisExp.addData('welcome.start', globalClock.getTime(format='float'))

# Draw welcome screen with Pokémon images
welcome_text = visual.TextStim(win=win, text="Welcome to the Pokémon Party game!", font='Arial', units='deg', pos=(0, 0), height=1.2, wrapWidth=1700, 
    color='black', colorSpace=CLR_SPC)
welcome_pokemon = {
    "Bulbasaur":  ((10, -5), (5, 5)),
    "Pikachu":    ((5, -5),  (5, 5)),
    "Squirtle":   ((0, -5),  (5, 5)),
    "Charmander": ((-5, -5), (5, 5)),
    "Krabby":     ((-10, -5), (5, 5)),
    "Magikarp":   ((-10, 5), (5, 5)),
    "Pidgey":     ((-5, 5),  (5, 5)),
    "Metapod":    ((0, 5),   (5, 5)),
    "Eevee":      ((5, 5),   (5, 5)),
    "Raticate":   ((10, 5),  (5, 5))} 

welcome_text.draw()
for name, (pos, size) in welcome_pokemon.items():
    stim = pokemon_dict[name]
    stim.pos = pos
    stim.size = size
    stim.draw()
win.flip()

# Reset size and position of all pokemon 
for pokemon in pokemon_dict:
    pokemon_dict[pokemon].size = POKEMON_SIZE
    pokemon_dict[pokemon].pos = (0,0)

# Wait for space or escape and clear the window
keys = event.waitKeys(keyList=['space', 'escape'])
if 'escape' in keys:
    terminate_task()
win.flip() 

###### EXPERIMENT #####################################################################################################################

# Generate or load the RSVPs, SIM onsets, and and peripheral stimulus grids for all runs
exp_rsvps = generate_rsvps()
exp_sim_onsets = generate_sim_onsets()
exp_trial_grids = assign_grids(feat_cond)

# Run requested FIX runs (1, 2, or 3)
if any(r in run_list for r in [1, 2, 3]):
    attention_cond = 'FIX'
    show_instructions(feat_cond, attention_cond)
    ready_screen()
    for run in run_list:
        if run <= RUNS_PER_COND:
            experiment_run(feat_cond, attention_cond, run, exp_rsvps[run], exp_sim_onsets[run], exp_trial_grids[run])

# Run requested COV runs (4, 5, or 6)
if any(r in run_list for r in [4, 5, 6]):
    attention_cond = 'COV'
    show_instructions(feat_cond, attention_cond)
    ready_screen()
    for run in run_list:
        if run >= RUNS_PER_COND:
            experiment_run(feat_cond, attention_cond, run, exp_rsvps[run], exp_sim_onsets[run], exp_trial_grids[run])

###### END EXPERIMENT ##################################################################################################################

# Draw thank you text 
thanks_text = visual.TextStim(win=win, text="Thanks for coming to the Pokémon Party!", font='Arial', units='deg', pos=(0, 0), height=1.2, wrapWidth=1700, 
    color='black', colorSpace=CLR_SPC)
thanks_text.draw()
win.flip()

# Wait until space is pressed to terminate task (saves date and closes window)
keys = event.waitKeys(keyList=['space'])
terminate_task() 