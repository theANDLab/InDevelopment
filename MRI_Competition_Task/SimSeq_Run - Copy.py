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
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
from psychopy.tools.monitorunittools import pix2deg
logging.console.setLevel(logging.ERROR)

###### PARAMETERS ######################################################################################################################

EYETRACKER_OFF = True # Set to True to run the script without eyetracking

# Initialize the global clock and keyboard
globalClock = core.Clock()
kb = keyboard.Keyboard(clock = globalClock)

# Practice 
PRACTICE_RUNS = 1 # a practice block consists of: blank, sim, seq, blank 
PRAC_SEQ_BLOCKS = 1 # how many seq blocks per practice run
PRAC_SIM_BLOCKS = 1 # how many sim blocks per practice run
PRAC_RSVP_RATE = 0.5 #sec; duration of RSVP pokemon presentation (slower than exp)
PRAC_BLOCK_DESIGN = [('RVF','SIM'),('LVF','SEQ')]

# Experiment
ATTENTION_CONDS = ['FIX', 'COV']
BLOCK_DESIGN = [('RVF','SIM'),('LVF','SEQ'),('RVF','SEQ'),('LVF','SIM'),('RVF','SEQ'),('LVF','SIM'),
                ('RVF','SIM'),('LVF','SEQ'),('RVF','SIM'),('LVF','SEQ'),('RVF','SEQ'),('LVF','SIM')]
                
NUM_RUNS = 6 # how many exp runs in one feature condition
RUNS_PER_COND = int(NUM_RUNS//len(ATTENTION_CONDS)) # equal number of FIX and COV runs per feature condition (runs 1,2,3 are FIX; 4,5,6 are COV)
NUM_SIM_BLOCKS = 6 # per run
NUM_SEQ_BLOCKS = 6 # per run
NUM_BLANK_BLOCKS = 2 # one before and one after each run (including practice run)
NUM_TRIALS = 3 # per block

# Timing
BLANK_BLOCK_DURATION = 16 # seconds
PERIPH_STIM_DURATION = 1 # seconds
RSVP_RATE = 0.25 #sec; duration of RSVP pokemon presentation
TRIAL_DURATION = PERIPH_STIM_DURATION*4 # sec
PSTIM_TARGET_FREQ = [1,3] # pstim color targets will occur every 1-3 trials
POKEMON_TARGET_FREQ = [15,30] # pokemon targets will occur every 15-30 pokemon (3.75-7.5s) in the RSVP
RESPONSE_WINDOW = 1.5 # sec; responses after this will be coded as FAs

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
SCANNER_KEY = '='

###### SETUP ###########################################################################################################################

pokemon_names = ["Bulbasaur", "Pikachu", "Squirtle", "Charmander", "Magikarp", "Raticate", "Pidgey",
    "Metapod", "Jigglypuff", "Butterfree", "Psyduck", "Caterpie", "Krabby",
    "Haunter", "Vulpix", "Eevee", "Sandshrew", "Wartortle", "Charmeleon", "Clefairy",
    "Ponyta", "Mankey"]
    
feat_conds = ['color', 'motion', 'color-motion'] # which feature conditions to display
    
exp_name = 'SimSeq'
exp_info = {
    'Participant ID': '9999', 
    'Session': '001',
    'Pokemon': pokemon_names,
    'Condition': feat_conds,
    'Runs': 'FIX: 1,2,3; COV: 4,5,6'
}
while True:
    dlg = gui.DlgFromDict(dictionary=exp_info, title=exp_name)
    if dlg.OK == False:
        core.quit()
        sys.exit()

    # write edf filename
    edf_filename = f"{exp_info['Participant ID']}"

    # check if the filename is valid
    allowed_char = ascii_letters + digits + '_'
    if not all([c in allowed_char for c in edf_filename]):
        raise ValueError('ERROR: Invalid EDF filename. Enter only letters, digits, or underscores.')
    elif len(edf_filename) > 8:
        raise ValueError("ERROR: Invalid EDF filename: participant ID must be ≤5 characters.")
    else:
        break

# Get and print task variabls
target_pokemon = exp_info['Pokemon'].strip().capitalize() 
feat_cond = exp_info['Condition']
run_list = [int(x.strip()) for x in exp_info['Runs'].split(',')] # which runs to display
sub_id = exp_info['Participant ID']

print('###############################################################################\n')
print(f'Performing runs {run_list} of {feat_cond} condition.')
print('Target pokemon:', target_pokemon)
print('Target color:', TARGET_COLOR)
print('\n###############################################################################')

# Establish data output directory and output file column order
time_str = time.strftime("%m_%d_%Y", time.localtime())
root_dir = os.path.dirname(os.path.abspath(__file__))
participant_folder = os.path.join(root_dir, 'data', f"{exp_info['Participant ID']}_{exp_name}_Session{exp_info['Session']}_{time_str}")
os.makedirs(participant_folder, exist_ok=True)
edf_path = os.path.join(participant_folder, f"{edf_filename}.EDF") # file for eyetracker data

trials_filename = os.path.join(participant_folder, f"trials_{exp_info['Participant ID']}_Session{exp_info['Session']}")
rsvp_filename = os.path.join(participant_folder, f"rsvp_{exp_info['Participant ID']}_Session{exp_info['Session']}")

# Create an experiment handler to manage the data file and set the column order
thisExp = data.ExperimentHandler(name=exp_name, version='', extraInfo=exp_info,
                                runtimeInfo=None, originPath=os.path.abspath(__file__),
                                savePickle=True, saveWideText=True,
                                dataFileName=trials_filename)

rsvpExp = data.ExperimentHandler(name='rsvp',extraInfo=exp_info,
                                savePickle=True, saveWideText=True,
                                dataFileName=rsvp_filename)
                                
column_order = ['welcome.start','instructions.start', 'instructions.end', 'prac_blank_block.start', 'prac_blank_block.end','blank_block.start', 'blank_block.end', 'feat_cond', 'run', 'attention_cond','block',
    'presentation_cond', 'vf', 'trial', 'rsvp_seq', 'pstim_colors', 'pstim_angles','trial.start', 'pstim.onset', 
    'target_shown', 'target.onset', 'press_times', 'rts', 'keypresses', 'hit']
    
rsvp_column_order = ['blank_block.start', 'rsvp_seq', 'run', 'block', 'trial']

for col in column_order:
    thisExp.addData(col, '')
    
for col in rsvp_column_order:
    rsvpExp.addData(col, '')
    
# Mark the experiment as started
exp_info['expDate'] = data.getDateStr(format = '%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6)
thisExp.status = STARTED

# Window setup (will need to be adjusted to match the MRI monitor)
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

# Get monitor's refresh rate
frame_rate = win.getActualFrameRate()
if frame_rate is not None:
    frame_dur = 1.0 / frame_rate
exp_info['frameRate'] = frame_rate

# Load Pokémon images from folder
image_dir = os.path.join(root_dir, 'images')
pokemon_dict = {name: visual.ImageStim(win, name=name, image=os.path.join(image_dir, f"{i+1:03}.png"))
    for i, name in enumerate(pokemon_names)} # dictionary of the pokemon where the key is their name and the values are the ImageStim object

# Calculate the coordinates of the center of the peripheral grid based on the inputted parameters
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

# -------------- Connect to the EyeLink Host PC; based on sample scripts from SR Research
if EYETRACKER_OFF:
    el_tracker = pylink.EyeLink(None) # Need this line so code will run when not using eyetracking
else:
    try:
        el_tracker = pylink.EyeLink("100.1.1.1")
    except RuntimeError as error:
        print('ERROR:', error)
        core.quit()
        sys.exit()
        
# Open the EDF data file on the Host PC
edf_file = edf_filename + ".EDF"
try:
    el_tracker.openDataFile(edf_file)
except RuntimeError as err:
    print('ERROR:', err)
    # close the link if we have one open
    if el_tracker.isConnected():
        el_tracker.close()
    core.quit()
    sys.exit()

# Add header text to EDF file for data viewing
preamble_text = 'RECORDED BY %s' % os.path.basename(__file__)
el_tracker.sendCommand("add_file_preamble_text '%s'" % preamble_text)

# Put the tracker in offline mode before we change tracking parameters
el_tracker.setOfflineMode()

# Get eyetracker version/model; EyeLink 1000 Plus is version 5
eyelink_ver = 0  # Set to 0 so code will run when not using eyetracking
if not EYETRACKER_OFF:
    vstr = el_tracker.getTrackerVersionString()
    eyelink_ver = int(vstr.split()[-1].split('.')[0])

# Set what eye events to save in the EDF file and make available over the link, include everything by default
file_event_flags = 'LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT'
link_event_flags = 'LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON,FIXUPDATE,INPUT'
# Set what sample data to save in the EDF data file and to make available over the link
if eyelink_ver > 3:
    file_sample_flags = 'LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT'
    link_sample_flags = 'LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT'
else: # For when running without eyetracking
    file_sample_flags = 'LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT'
    link_sample_flags = 'LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT'
el_tracker.sendCommand("file_event_filter = %s" % file_event_flags)
el_tracker.sendCommand("file_sample_data = %s" % file_sample_flags)
el_tracker.sendCommand("link_event_filter = %s" % link_event_flags)
el_tracker.sendCommand("link_sample_data = %s" % link_sample_flags)

# Choose a calibration type (HV = horizontal/vertical) number is how many points on the screen
el_tracker.sendCommand("calibration_type = HV5")

# Get the screen resolution used by PsychoPy
scn_width, scn_height = win.size # in retina pixels

# Pass the display pixel coordinates (left, top, right, bottom) to the tracker
el_coords = "screen_pixel_coords = 0 0 %d %d" % (scn_width - 1, scn_height - 1)
el_tracker.sendCommand(el_coords)

# Write a DISPLAY_COORDS message to the EDF file
# Data Viewer needs this piece of info for proper visualization
dv_coords = "DISPLAY_COORDS  0 0 %d %d" % (scn_width - 1, scn_height - 1)
el_tracker.sendMessage(dv_coords)

# Configure a graphics environment (genv) for tracker calibration
genv = EyeLinkCoreGraphicsPsychoPy(el_tracker, win)

# Set visuals for calibration routine
foreground_color = (-1, -1, -1)
background_color = win.color
genv.setCalibrationColors(foreground_color, background_color)
genv.setTargetType('picture')
genv.setPictureTarget(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'andy_fixation.png')) 

# Beeps to play during calibration, validation and drift correction
# parameters: target, good, error
#     target -- sound to play when target moves
#     good -- sound to play on successful operation
#     error -- sound to play on failure or interruption
# Each parameter could be ''--default sound, 'off'--no sound, or a wav file
genv.setCalibrationSounds('boing.wav', '', '')

# Request Pylink to use the PsychoPy window we opened above for calibration
pylink.openGraphicsEx(genv)
logging.info(f"Graphics environment set up: {genv}")

###### FUNCTIONS #######################################################################################################################

def terminate_task():
    """ Saves data and closes the window."""
    
    thisExp.nextEntry()
    thisExp.addData('experiment.end', globalClock.getTime(format='float'))
    thisExp.saveAsWideText(trials_filename + '.csv', delim='auto')
    thisExp.saveAsPickle(trials_filename)
    logging.flush()
    
    # Clear psychopy window
    if win is not None:
        win.clearAutoDraw()
        win.flip()
        
    # Mark experiment as finished
    thisExp.status = FINISHED
    print("Experiment ended.")
    thisExp.abort()
        
    if not EYETRACKER_OFF:
        # Disconnect eyetracker
        el_tracker = pylink.getEYELINK()
        if el_tracker.isConnected():
            # Put tracker in Offline mode
            el_tracker.setOfflineMode()
            # Clear the Host PC screen and wait for 500 ms
            el_tracker.sendCommand('clear_screen 0')
            pylink.msecDelay(500)
            # Close the edf data file on the Host
            el_tracker.closeDataFile()
            # Print a file transfer message
            print('EDF data is transferring from EyeLink Host PC...')
            # Download the EDF data file from the Host PC to a local data folder
            try:
                el_tracker.receiveDataFile(edf_file, edf_path)
                print(f"EDF file saved to: {edf_path}")
            except RuntimeError as error:
                print('ERROR downloading EDF file:', error)
            # Close the link to the tracker
            el_tracker.close()

    print("Task ended.")
    thisExp.abort() # or data files will save again on exit
    win.close()
    core.quit()
    sys.exit()
    
def drift_check():
    """ Performs a drift check. Allows for recalibration during the exp. """
    
    # the doDriftCorrect() function requires target position in integers
    # the last two arguments:
    # draw_target (1-default, 0-draw the target then call doDriftCorrect)
    # allow_setup (1-press ESCAPE to recalibrate, 0-not allowed)
    while not EYETRACKER_OFF:
        # terminate the task if no longer connected to the tracker
        if (not el_tracker.isConnected()) or el_tracker.breakPressed():
            terminate_task()
            
        # drift-check and re-do camera setup if ESCAPE is pressed
        try:
            error = el_tracker.doDriftCorrect(int(scn_width/2.0),int(scn_height/2.0), 1, 1)
            # break following a success drift-check
            if error is not pylink.ESC_KEY:
                break
        except:
            pass
    
def get_eye_used(el_tracker):
    """ Gets eye used. Returns 0 for left, 1 for right, None if eye data cannot be collected.t"""
    if EYETRACKER_OFF:
        return 0 # For when running without eyetracking
        
    if el_tracker is not None and el_tracker.isConnected():
        eye = el_tracker.eyeAvailable()
        if eye in [0, 1]:
            el_tracker.sendMessage(f"EYE_USED {eye} {'RIGHT' if eye == 1 else 'LEFT'}")
            return eye
        elif eye == 2: # binocular vision defaults to left eye
            el_tracker.sendMessage("EYE_USED 0 LEFT")
            return 0
            
    print("ERROR: EyeLink not connected or invalid eye")
    return None
    
def is_gaze_within_bounds(el_tracker, eye_used):
    """ Check that gaze is within GAZE_BOUND, return False if not. """
    while True:
        sample = el_tracker.getNewestSample()
        
        if eye_used == 0 and sample.isLeftSample():
            eye_data = sample.getLeftEye()
        elif eye_used == 1 and sample.isRightSample():
            eye_data = sample.getRightEye()
        
        gaze = eye_data.getGaze()
        pupil = eye_data.getPupilSize()
        
        dx = gaze[0] - scn_width/2
        dy = scn_height/2 - gaze[1]
        pixels = np.array([dx, dy])
        dx_deg, dy_deg = pix2deg(pixels, monitor=Eizo)
        
        return abs(dx_deg) <= GAZE_BOUND and abs(dy_deg) <= GAZE_BOUND
    
def generate_blank_rsvps():
    """" 
    Generates lists of pokemon to be presented during all blank blocks in the experiment. 
    
    Returns two dictionaries, one for practice runs and one for experiment runs. 
    Keys are run numbers. Values for each key are 2 lists of pokemon names that serve as the RSVP stimulus 
    presentation lists for the blank blocks of that run (one before and one after). 
    
    Checks if dictionaries have already been generated for this subID. 
    """
    
    # Check if dictionaries exist in pt folder. If so, load and return it
    filename = os.path.join(participant_folder, f"{sub_id}_blank_rsvps.json")
    
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        prac_blank_rsvps = {int(k): v for k, v in data['prac_blank_rsvps'].items()}
        exp_blank_rsvps = {int(k): v for k, v in data['exp_blank_rsvps'].items()}
        print(f"Loaded existing blanks RSVPs from {filename}")
        return prac_blank_rsvps, exp_blank_rsvps
    
    prac_blank_sequences = []
    prac_blanks = NUM_BLANK_BLOCKS * PRACTICE_RUNS # total number of blank blocks across all practice runs
    pokemon_per_prac = int(BLANK_BLOCK_DURATION // PRAC_RSVP_RATE)

    exp_blank_sequences = []
    exp_blanks = NUM_BLANK_BLOCKS * RUNS_PER_COND # total number of blank blocks across all experiment runs
    pokemon_per_blank = int(BLANK_BLOCK_DURATION // RSVP_RATE)
    
    # Create all of the lists of pokemon names for practice blank blocks. 
    for blank_block in range(prac_blanks):
        # Get times for target occurrences in this block
        target_indices = []
        next_target_idx = random.randint(*POKEMON_TARGET_FREQ)
        while next_target_idx < pokemon_per_prac:
            target_indices.append(next_target_idx)
            next_target_idx += random.randint(*POKEMON_TARGET_FREQ)
        # Build sequence, prevent back-to-back repeats
        sequence = []
        for idx in range(pokemon_per_prac):
            if idx in target_indices:
                sequence.append(target_pokemon)
            else:
                distractor_options = [p for p in pokemon_names if p != target_pokemon]
                if idx > 0:
                    distractor_options = [p for p in distractor_options if p != sequence[-1]] #ensure no back to back pokemon
                sequence.append(random.choice(distractor_options))
        prac_blank_sequences.append(sequence)
    
    # Create all of the lists of pokemon names for experiment blank blocks. 
    for blank_block in range(exp_blanks):
        # Get times for target occurrences in this block
        target_indices = []
        next_target_idx = random.randint(*POKEMON_TARGET_FREQ)
        while next_target_idx < pokemon_per_blank:
            target_indices.append(next_target_idx)
            next_target_idx += random.randint(*POKEMON_TARGET_FREQ)
        # Build sequence, prevent back-to-back repeats
        sequence = []
        for idx in range(pokemon_per_blank):
            if idx in target_indices:
                sequence.append(target_pokemon)
            else:
                distractor_options = [p for p in pokemon_names if p != target_pokemon]
                if idx > 0:
                    distractor_options = [p for p in distractor_options if p != sequence[-1]] #ensure no back to back pokemon
                sequence.append(random.choice(distractor_options))
        exp_blank_sequences.append(sequence)

    # Assign 2 rsvp sequences to each stimulus set number
    prac_blank_rsvps = {}
    exp_blank_rsvps = {}
    idx = 0
    for run in range(PRACTICE_RUNS):
        prac_blank_rsvps[run+1] = [prac_blank_sequences[idx], prac_blank_sequences[idx + 1]]
        idx +=NUM_BLANK_BLOCKS
    idx = 0
    for run in range(RUNS_PER_COND):
        exp_blank_rsvps[run+1] = [exp_blank_sequences[idx], exp_blank_sequences[idx + 1]]
        idx +=NUM_BLANK_BLOCKS
        
    # Blank blocks for runs 4,5,6 are the same as for 1,2,3
    exp_blank_rsvps[4] = exp_blank_rsvps[1]
    exp_blank_rsvps[5] = exp_blank_rsvps[2]
    exp_blank_rsvps[6] = exp_blank_rsvps[3]
    
    # Save dictionaries to file for future loading
    data = {
        'prac_blank_rsvps': {str(k): v for k, v in prac_blank_rsvps.items()},
        'exp_blank_rsvps': {str(k): v for k, v in exp_blank_rsvps.items()}
    }
    with open(filename, 'w') as f:
        json.dump(data, f)
    print(f"Generated and saved blank RSVPs to {filename}")
    
    return prac_blank_rsvps, exp_blank_rsvps
    
def generate_trial_rsvps():
    """ 
    Generates pokemon RSVP lists for all trials in the experiment.
    
    Returns two dictionaries, one for practice runs and one for experiment runs. 
    Keys are run numbers with x number of values. X is the total number of trials in that run. 
    Values are lists of pokemon names that serve as the RSVP stimulus lists for that trial. 
    
    Example: prac_trial_rsvps[1][0] is the list of pokemon to be presented in the
    first trial of the first practice run. Trials are 0-indexed.
    
    Checks if dictionaries have already been generated for this subID.
    """
    
    # Check if dictionaries exist in pt folder. If so, load and return it
    filename = os.path.join(participant_folder, f"{sub_id}_trial_rsvps.json")
    
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        prac_trial_rsvps = {int(k): v for k, v in data['prac_trial_rsvps'].items()}
        exp_trial_rsvps = {int(k): v for k, v in data['exp_trial_rsvps'].items()}
        print(f"Loaded existing trial RSVPs from {filename}")
        return prac_trial_rsvps, exp_trial_rsvps
        
    exp_seq_list = []
    pokemon_per_trial = int(TRIAL_DURATION // RSVP_RATE)
    total_pokemon = int((NUM_TRIALS*(NUM_SIM_BLOCKS+NUM_SEQ_BLOCKS))*pokemon_per_trial) # in one exp run
    
    prac_seq_list = []
    pokemon_per_prac = int(TRIAL_DURATION // PRAC_RSVP_RATE)
    practice_pokemon = int(NUM_TRIALS*(PRAC_SEQ_BLOCKS+PRAC_SIM_BLOCKS)*pokemon_per_prac)
    
    for run in range(PRACTICE_RUNS):
        target_pokemon_idx = []
        next_target_idx = random.randint(*POKEMON_TARGET_FREQ)
        while next_target_idx < practice_pokemon:
            target_pokemon_idx.append(next_target_idx)
            next_target_idx += random.randint(*POKEMON_TARGET_FREQ)
        run_sequence = []
        for rsvp_idx in range(practice_pokemon):
            if rsvp_idx in target_pokemon_idx:
                run_sequence.append(target_pokemon)
            else:
                distractors = [p for p in pokemon_names if p!= target_pokemon]
                if rsvp_idx > 0:
                    distractors = [p for p in distractors if p != run_sequence[-1]]
                run_sequence.append(random.choice(distractors))
        prac_seq_list.append(run_sequence)
    
    prac_trial_rsvps = {}
    for run, run_sequence in enumerate(prac_seq_list):
        trial_sequences = [run_sequence[i:i + pokemon_per_prac] for i in range(0, len(run_sequence), pokemon_per_prac)]
        prac_trial_rsvps[run+1] = trial_sequences
            
    for run in range(RUNS_PER_COND):
        target_pokemon_idx = []
        next_target_idx = random.randint(*POKEMON_TARGET_FREQ)
        while next_target_idx < total_pokemon:
            target_pokemon_idx.append(next_target_idx)
            next_target_idx += random.randint(*POKEMON_TARGET_FREQ) # list of all target indices for the run
        # Create the sequence for the entire run
        run_sequence = []
        for rsvp_idx in range(total_pokemon):
            if rsvp_idx in target_pokemon_idx:
                run_sequence.append(target_pokemon)
            else:
                distractors = [p for p in pokemon_names if p!= target_pokemon]
                if rsvp_idx > 0:
                    distractors = [p for p in distractors if p != run_sequence[-1]]
                run_sequence.append(random.choice(distractors))
        exp_seq_list.append(run_sequence)

    exp_trial_rsvps = {}
    for run, run_sequence in enumerate(exp_seq_list):
        trial_sequences = [run_sequence[i:i + pokemon_per_trial] for i in range(0, len(run_sequence), pokemon_per_trial)]
        exp_trial_rsvps[run+1] = trial_sequences
        
    # Trial rsvps for runs 1,2,3 are the same for 4,5,6
    exp_trial_rsvps[4] = exp_trial_rsvps[1]
    exp_trial_rsvps[5] = exp_trial_rsvps[2]
    exp_trial_rsvps[6] = exp_trial_rsvps[3]
    
    # Save dictionaries for future reloading
    data = {
        'prac_trial_rsvps': {str(k): v for k, v in prac_trial_rsvps.items()},
        'exp_trial_rsvps': {str(k): v for k, v in exp_trial_rsvps.items()}
    }
    with open(filename, 'w') as f:
        json.dump(data,f)
    print(f"Generated and saved trial RSVPs to {filename}")
    
    return prac_trial_rsvps, exp_trial_rsvps
    
def generate_sim_onsets():
    """
    Generates pseudo-random onset times for peripheral stimulus grid presentation in SIM trials. 
    
    Returns two dictionaries; one for practice runs and one for exp runs. Each dictionary has
    x keys, one for each run.Values are lists of the onset times for all SIM trials of that run.
    
    Example: prac_sim_onsets[1] is list of integers, one for each SIM trial of the first practice run. 
    
    Onsets will never be back to back (i.e., if the grid was presented during the last second of the previous trial,
    it will not be presented during the first second of the current trial). 
    """
    # Check whether dictionaries have already been generated for this subID
    filename = os.path.join(participant_folder, f"{sub_id}_sim_onsets.json")
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        prac_sim_onsets = {int(k): v for k, v in data['prac_sim_onsets'].items()}
        exp_sim_onsets = {int(k): v for k, v in data['exp_sim_onsets'].items()}
        print(f"Loaded existing blank RSVPs from {filename}")
        return prac_sim_onsets, exp_sim_onsets
    
    prac_sim_onsets = {}
    for prac_run in range(1, PRACTICE_RUNS+1):
        trial_onsets = []
        for trial_idx in range(PRAC_SIM_BLOCKS*NUM_TRIALS):
            if trial_idx != 0 and trial_onsets[trial_idx-1] == 3:
                trial_onsets.append(random.choice([1, 2, 3]))
            else:
                trial_onsets.append(random.choice([0, 1, 2, 3]))
        prac_sim_onsets[prac_run] = trial_onsets
    
    exp_sim_onsets = {}
    for run in range(1, RUNS_PER_COND+1):
        trial_onsets = []
        for trial_idx in range(NUM_SIM_BLOCKS*NUM_TRIALS):
            if trial_idx != 0 and trial_onsets[trial_idx-1] == 3:
                trial_onsets.append(random.choice([1, 2, 3]))
            else:
                trial_onsets.append(random.choice([0, 1, 2, 3]))
        exp_sim_onsets[run] = trial_onsets
        
    # SIM onsets for runs 1,2,3 are the same for 4,5,6
    exp_sim_onsets[4] = exp_sim_onsets[1]
    exp_sim_onsets[5] = exp_sim_onsets[2]
    exp_sim_onsets[6] = exp_sim_onsets[3]
    
    # Save SIM onsets to file for future reloading
    data = {
        'prac_sim_onsets': {str(k): v for k, v in prac_sim_onsets.items()},
        'exp_sim_onsets': {str(k): v for k, v in exp_sim_onsets.items()}
    }
    with open(filename, 'w') as f:
        json.dump(data, f)
    print(f"Generated and saved blank RSVPs to {filename}")
    
    return prac_sim_onsets, exp_sim_onsets
    
def assign_grids(feat_cond):
    """ 
    Returns two dictionaries (keys: run) of lists of each trial's peripheral stim parameters 
    (colors, angles). One for practice runs, one for experiment runs.
        
    Example:
        Calling exp_trial_grids[1][32] will give the grid layout (list of colors and angles) for 
        the 33rd trial in the first experiment run. 
    """
    # Check whether dictionaries have already been generated for this subID
    filename = os.path.join(participant_folder, f"{sub_id}_grids.json")
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        prac_trial_grids = {int(k): v for k, v in data['prac_trial_grids'].items()}
        exp_trial_grids = {int(k): v for k, v in data['exp_trial_grids'].items()}
        print(f"Loaded existing grids from {filename}")
        return prac_trial_grids, exp_trial_grids

    trials_per_run = NUM_TRIALS * (NUM_SEQ_BLOCKS + NUM_SIM_BLOCKS)
    prac_trials = NUM_TRIALS * (PRAC_SEQ_BLOCKS + PRAC_SIM_BLOCKS)
    prac_trial_grids = {}
    exp_trial_grids = {}

    #-------Color feature condition
    all_color_combos = list(itertools.combinations(PERIPHERAL_STIM_COLORS.keys(), NUM_PSTIMS))
    all_color_configs = []
    for combo in all_color_combos:
        all_color_configs.extend(itertools.permutations(combo))
    random.shuffle(all_color_configs)
    
    rvf_target_colors = [grid for grid in all_color_configs if TARGET_COLOR in grid and grid[2] == TARGET_COLOR]
    lvf_target_colors = [grid for grid in all_color_configs if TARGET_COLOR in grid and grid[3] == TARGET_COLOR]
    colors_without_target = [grid for grid in all_color_configs if TARGET_COLOR not in grid]
    
    #-------Motion and color-motion condition
    all_angle_combos = list(itertools.combinations(ANGLES, NUM_PSTIMS))
    all_angle_configs = []
    for combo in all_angle_combos:
        all_angle_configs.extend(itertools.permutations(combo))
    random.shuffle(all_angle_configs)
    
    rvf_target_angles = [c for c in all_angle_configs if c[2] == TARGET_ANGLE]
    lvf_target_angles = [c for c in all_angle_configs if c[3] == TARGET_ANGLE]
    angles_without_target = [c for c in all_angle_configs if TARGET_ANGLE not in c]

    #-------Helper function to assign grids for a single run
    def assign_run_grids(num_trials, block_design):
        target_trials = []
        target_trial_idx = random.randint(*PSTIM_TARGET_FREQ)
        while target_trial_idx <= num_trials:
            target_trials.append(target_trial_idx)
            target_trial_idx += random.randint(*PSTIM_TARGET_FREQ)

        run_assignments = []
        used_color_configs = set()
        used_angle_configs = set()

        for trial_idx in range(num_trials):
            block_num = (trial_idx // NUM_TRIALS) % len(block_design)
            block_vf = block_design[block_num][0]
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
                    color_pool = colors_without_target + (rvf_target_colors if block_vf == 'RVF' else lvf_target_colors)
                    color_available = [c for c in color_pool if c not in used_color_configs] or colors_without_target
                    angle_available = [a for a in angles_without_target if a not in used_angle_configs] or angles_without_target
                chosen_colors = random.choice(color_available)
                chosen_angles = random.choice(angle_available)
                used_color_configs.add(chosen_colors)
                used_angle_configs.add(chosen_angles)
                run_assignments.append({'colors': chosen_colors, 'angles': chosen_angles})

        return run_assignments

    #-------Assign grids for practice runs
    for run in range(1, PRACTICE_RUNS + 1):
        prac_trial_grids[run] = assign_run_grids(prac_trials, PRAC_BLOCK_DESIGN)

    #-------Assign grids for experiment runs
    for run in range(1, RUNS_PER_COND + 1):
        exp_trial_grids[run] = assign_run_grids(trials_per_run, BLOCK_DESIGN)

    # Runs 4,5,6 reuse the same grids as runs 1,2,3
    exp_trial_grids[4] = exp_trial_grids[1]
    exp_trial_grids[5] = exp_trial_grids[2]
    exp_trial_grids[6] = exp_trial_grids[3]

    # Save to file for future reloading
    data = {
        'prac_trial_grids': {str(k): v for k, v in prac_trial_grids.items()},
        'exp_trial_grids': {str(k): v for k, v in exp_trial_grids.items()}
    }
    with open(filename, 'w') as f:
        json.dump(data, f)
    print(f"Generated and saved grids to {filename}")

    return prac_trial_grids, exp_trial_grids

    
def create_trial_dicts(trial_grids, trial_rsvps, run_sim_onsets, practice = False):
    """
    Returns a list of trial dictionaries for all of the trials in a run. Call at the start of each run. 

    Args:
        trial_grids (list): list of dictionaries, one for each trial's grid parameters
        trial_rsvps (list): rsvp sequences for all trials in this run
        
    Returns:
        trial_dicts (list): Each value is a trial dictionary including:
            'block_num', 'trial_num', 'visual_field', 'presentation_cond', 
            'peripheral_grid', 'rsvp_seq'
            
    Example:
        Calling create_trial_dicts(trial_grids, trial_rsvps, run_sim_onsets)[32] returns the trial dictionary of
        the 33rd trial in the run
    """    
    sim_onsets = deque(run_sim_onsets) # .popleft is much faster with deque object than list
    trial_dicts = []
    if practice:
        design = PRAC_BLOCK_DESIGN
    else:
        design = BLOCK_DESIGN
        
    for block_idx, (visual_field, present_cond) in enumerate(design):
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
    # Show instructions and target for the attention condition
    thisExp.addData('instructions.start', globalClock.getTime(format='float'))
    
    # Instructions based on target in cov condition
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

    # Fix instructions are the same throughout all feature conditions
    fix_instructions_text = visual.TextStim(win=win, text=(
        "There's a Pokémon Party happening right now, and the Pokémon are playing hide and seek!\n\n"
        f"The Pokémon are trying to find {target_pokemon}! Can you help them? {target_pokemon} will show up like this:\n\n\n"
        f"\n\nPress the button as fast as you can every time you see {target_pokemon}.\n\n\n"
        "Ready to start playing?"), font='Arial', units='deg', pos=(0, 0), height=1, wrapWidth=1400, 
        color='black', colorSpace=CLR_SPC)
    
    # Draw instructions components on the screen
    if attention_cond == "FIX":
        fix_instructions_text.draw()
        pokemon_dict[target_pokemon].pos = (0, 0) 
        pokemon_dict[target_pokemon].size = [1.5, 1.5]
        pokemon_dict[target_pokemon].draw()
        win.flip()
        keys = event.waitKeys(keyList=['space', 'escape'])
        
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
    
def show_feedback(feat_cond, attention_cond):
    """ Function to display feedback based on feat and attention conds."""
    
    thisExp.addData('feedback.start', globalClock.getTime(format='float'))
    
    end_text = visual.TextStim(win=win, text=(), font='Arial', units='deg', pos=(0, 0), height=1.2, wrapWidth=1700, 
        color='black', colorSpace=CLR_SPC)
        
    if attention_cond == 'FIX':
        end_text.text = (f"Great job finding {target_pokemon}!")
        pokemon_dict[target_pokemon].pos = (0, -5) 
        pokemon_dict[target_pokemon].size = (5,5)
        pokemon_dict[target_pokemon].draw()
        win.flip()
        keys = event.waitKeys(keyList=['space', 'escape'])
        
        # Reset size and position of the target pokemon
        pokemon_dict[target_pokemon].size = POKEMON_SIZE
        pokemon_dict[target_pokemon].pos = (0,0)
        
        if 'escape' in keys:
            terminate_task()
            
    elif attention_cond == 'COV':
        if feat_cond == 'motion':
            target_pstim = visual.Circle(win, name = TARGET_COLOR, pos = (0,-5), radius = 1 , units = 'deg', 
                anchor = 'center', fillColor='black', lineColor='black')
        elif feat_cond == 'color-motion':
            target_pstim = visual.Circle(win, name = TARGET_COLOR, pos = (0,-5), radius = 1 , units = 'deg', anchor = 'center', 
                fillColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], lineColor=PERIPHERAL_STIM_COLORS[TARGET_COLOR], colorSpace = CLR_SPC)
        end_text.text = (f"Great job feeding the Pokémon the circles!")
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
                
    thisExp.addData('feedback.start', globalClock.getTime(format='float'))
    
def ready_screen():
    """ Text screen before the start of the experiment trials."""
    
    thisExp.addData('ready_screen.start', globalClock.getTime(format='float'))

    ready_text = visual.TextStim(win=win, text=("Ready to start the real game?"), font='Arial', units='deg', pos=(0, 0), height=1.2, wrapWidth=1700, 
        color='black', colorSpace=CLR_SPC)
        
    ready_text.draw()
    win.flip()
    
    keys = event.waitKeys(keyList=['space', 'escape'])
    if 'escape' in keys:
        terminate_task()
        
    thisExp.addData('ready_screen.end', globalClock.getTime(format='float'))

def run_blank_block(rsvp, run_num, blank_num, attn_cond, feat_cond, last_target_onset, practice=False): 
    """ Function to run a single blank block. """
    
    # Add data to data file
    block_start = globalClock.getTime(format='float')
    rsvpExp.addData('blank_block.start', block_start)
    rsvpExp.addData('rsvp_seq', rsvp)
    
    # -------------------- Eyetracker Setup ---------------------------------
    # Esure tracker is ready to receive commands
    el_tracker = pylink.getEYELINK()
    el_tracker.setOfflineMode()
    el_tracker.sendCommand('clear_screen 0')
    
    # Print trial number on eyelink host monitor and output console
    status_msg = 'Blank block'
    el_tracker.sendMessage('TRIALID BLANK')

    # Send status message to host PC
    el_tracker.sendCommand("record_status_message '%s'" % status_msg)

    # put tracker in idle/offline mode before recording
    el_tracker.setOfflineMode()
    
    # Start recording
    try:
        el_tracker.startRecording(1, 1, 1, 1) # arguments: sample_to_file, events_to_file, sample_over_link, event_over_link (1-yes, 0-no)
    except RuntimeError as error:
        print("ERROR:", error)
    
    # Allocate time for the tracker to cache some samples
    pylink.pumpDelay(100) 
    
    # Get eye used
    eye_used = get_eye_used(el_tracker)
    if eye_used is None and not EYETRACKER_OFF:
        print(f"Could not get eye used on trial {trial_dict['trial_num']}.")
        
    # Mark trial if eyetracker disconnected
    error = el_tracker.isRecording()
    if error is not pylink.TRIAL_OK:
        el_tracker.sendMessage('tracker_disconnected')
        print("Tracker disconnected.")
        
    # ------------------------------------------------------------------------
    
    if practice:
        RATE = PRAC_RSVP_RATE
    else:
        RATE = RSVP_RATE
        
    next_pokemon_onset = block_start
    for current_pokemon in rsvp:
        rsvpExp.addData('feat_cond', feat_cond)
        if practice:
            rsvpExp.addData('run', f'prac{run_num}')
        else:
            rsvpExp.addData('run', f'exp{run_num}')
        rsvpExp.addData('block', 'blank')
        rsvpExp.addData('attention_cond', attn_cond)
        hit = 0
        is_target = (current_pokemon == target_pokemon)
        
        # Set pokemon location
        pokemon_dict[current_pokemon].pos = (0,0)

        # Reset onset variables
        pokemon_onset = None
        timer_end = next_pokemon_onset + RATE
        
        # while loop will draw pokemon and check for keypresses for RATE duration
        while globalClock.getTime() < timer_end:
            if not EYETRACKER_OFF:
                if not is_gaze_within_bounds(el_tracker, eye_used):
                    print(f"Fixation broken during blank block.")
                    rsvpExp.addData('fix.broken', 'fix.broken')
            
            # Add pokemon to drawing queue
            pokemon_dict[current_pokemon].draw()
            
            # Capture pokemon onset once during first window flip
            if pokemon_onset is None:
                win.callOnFlip(kb.clearEvents)
                def on_flip():
                    nonlocal pokemon_onset, last_target_onset
                    t = globalClock.getTime(format='float')
                    pokemon_onset = t
                    rsvpExp.addData('stim', current_pokemon)
                    rsvpExp.addData('stim.onset', t) # onset will be on window flip to capture exactly when pokemon was first displayed
                    if is_target:
                        last_target_onset = t
                        rsvpExp.addData('target.onset', t)
                win.callOnFlip(on_flip)
            win.flip()
            
            # Get all keys pressed during that 250ms window
            keys = kb.getKeys(keyList=[RESPONSE_KEY, 'escape'], waitRelease=False, clear=False)
            
            # Analyze key presses
            for key in keys:
                if key.name == 'escape':
                    terminate_task()
                    
                elif key.name == RESPONSE_KEY:
                    press_time = key.rt # relative to globalClock
                    rsvpExp.addData('press_time', press_time)
                    # Score against the most recent target onset
                    if last_target_onset is not None:
                        # compute rt since last target onset
                        rt = press_time - last_target_onset
                        rsvpExp.addData('rt', rt)
                        if rt <= RESPONSE_WINDOW:
                            hit = 1
                            rsvpExp.addData('hit', hit)
                            
        def offset_on_flip():
            t = globalClock.getTime(format='float')
            rsvpExp.addData('stim.offset', t)
            
        win.callOnFlip(offset_on_flip)
        win.flip(clearBuffer = True)
        rsvpExp.nextEntry()
        
        next_pokemon_onset += RATE
            
    # Save block data
    rsvpExp.addData('blank_block.end', globalClock.getTime(format='float'))
    rsvpExp.nextEntry()
    
    # Send trial data to EDF file
    
    el_tracker.sendMessage('!V CLEAR 128 128 128')
    
    # Stop recording between trials to decrease size of output file
    pylink.pumpDelay(100) # add 100 msec to catch final events before stopping
    el_tracker.stopRecording()
    
    # Send trial result message to mark the end of the trial
    el_tracker.sendMessage('TRIAL_RESULT %d' % pylink.TRIAL_OK)
    # Stop recording between trials to decrease size of output file
    pylink.pumpDelay(100) # add 100 msec to catch final events before stopping
    el_tracker.stopRecording()
    
    return last_target_onset

def run_trial(feat_cond, run, trial_dict, attention_cond, last_target_onset, practice = False):    
    """
    Function to run a single trial.

    Args:
       trial_dict(dict): Contains 'block_num', 'trial_num', 'visual_field', 'presentation_cond', 'peripheral_grid','rsvp_sequence'
       attention_cond (str): 'FIX' (looking for a pokemon) or 'COV' (looking for a colored circle)
       target_pokemon (str): Name of pokemon to look for
    """
    
    # Reset variables for trial
    kb.clearEvents()
    rsvp_idx = 0
    pstim_idx = 0
    hit = 0
    press_times = []
    rts = []
    target_onset_recorded = False
    
    # Non-slip timing using global clock
    trial_start = globalClock.getTime()
    trial_end = trial_start + TRIAL_DURATION
    
    # Get presentation condition
    is_seq_trial = trial_dict['presentation_cond'] == 'SEQ'
    is_sim_trial = trial_dict['presentation_cond'] == 'SIM'
    
    # Set up rsvp sequence
    rsvp_sequence = trial_dict['rsvp_seq']
    next_pokemon_onset = trial_start
    current_pokemon = None
    if practice:
        RATE = PRAC_RSVP_RATE
    else:
        RATE = RSVP_RATE

    # Extract peripheral stim locations
    vf = trial_dict['visual_field']
    if vf[0] == 'R':
        grid_positions = [rvf_topright, rvf_topleft, rvf_botleft, rvf_botright]
    elif vf[0] == 'L':
        grid_positions = [lvf_topright, lvf_topleft, lvf_botleft, lvf_botright]

    # Map peripheral stim colors, angles, and phases to their positions 
    pstim_to_draw = []
    pstim_info = {}
    pstim_colors = trial_dict['peripheral_grid']['colors']
    pstim_angles = trial_dict['peripheral_grid']['angles']
    phases = [0, 0, np.pi/2, np.pi/2] # half of the circles will start at middle of motion, and half will start at a peak
    random.shuffle(phases)
    for color, angle, pos, phase in zip(pstim_colors, pstim_angles, grid_positions, phases):
        if color == 'black':
            pstim = visual.Circle(win, name = f"{color}_{angle}", pos = pos, radius = PERIPHERAL_STIM_SIZE/2 , 
                units = 'deg', anchor = 'center', fillColor=color, lineColor=color, colorSpace=CLR_SPC)
        else:
            pstim = visual.Circle(win, name = f"{color}_{angle}", pos = pos, radius = PERIPHERAL_STIM_SIZE/2 , 
                units = 'deg', anchor = 'center', fillColor=PERIPHERAL_STIM_COLORS[color], lineColor=PERIPHERAL_STIM_COLORS[color], colorSpace=CLR_SPC)
        pstim_to_draw.append(pstim) # create a list of pstims to draw in this trial
        pstim_info[pstim.name] = {'angle': angle, 'base_pos': pstim.pos, 'phase': phase}
    pstim_onset_recorded_dict = {pstim.name: False for pstim in pstim_to_draw}
    pstim_offset_recorded_dict = {pstim.name: False for pstim in pstim_to_draw}
        
    # Save trial variables to the data file
    thisExp.addData('feat_cond', feat_cond)
    thisExp.addData('run', run)
    thisExp.addData('block', trial_dict['block_num'] )
    thisExp.addData('trial', trial_dict['trial_num'])
    thisExp.addData('attention_cond', attention_cond)
    thisExp.addData('presentation_cond', trial_dict['presentation_cond'])
    thisExp.addData('vf', vf[0])
    thisExp.addData('rsvp_seq', rsvp_sequence)
    thisExp.addData('pstim_colors', pstim_colors)
    thisExp.addData('pstim_angles', pstim_angles)
    thisExp.addData('pstim_phases', phases)
    thisExp.addData('trial.start', trial_start)

    # Set peripheral stim grid onsets
    if is_sim_trial:
        pstim_start_time = trial_dict['grid_onset'] # start time for entire grid in sim trials
        thisExp.addData('pstim.onset', pstim_start_time)
    else:
        pstim_onsets = [0, 1, 2, 3] # seconds; SEQ condition: each peripheral stimuli will be shown one at a time for 1s
        thisExp.addData('pstim.onset', 0) # first pstim is always at start of trial
    
    # Check whether target is present this trial
    if attention_cond == 'FIX': # target shown if target pokemon appears in rsvp during the trial
        target_shown = target_pokemon in rsvp_sequence 
    elif attention_cond == 'COV':
        if feat_cond == 'color': # target shown if target color appears in the grid this trial 
            target_shown = TARGET_COLOR in pstim_colors
        elif feat_cond == 'motion': # target shown if black circle closest to fix is moving vertically
            if vf[0] == 'L':
                target_shown = TARGET_ANGLE == pstim_angles[3]
            else:
                target_shown = TARGET_ANGLE == pstim_angles[2]
        elif feat_cond == 'color-motion': # target shown if circle closest to fix is moving vertically and is target color
            if vf[0] == 'L':
                target_shown = (TARGET_ANGLE == pstim_angles[3] and TARGET_COLOR == pstim_colors[3])
            else:
                target_shown = (TARGET_ANGLE == pstim_angles[2] and TARGET_COLOR == pstim_colors[2])
    
    # -------------------- Eyetracker Setup ---------------------------------
    # Esure tracker is ready to receive commands
    el_tracker = pylink.getEYELINK()
    el_tracker.setOfflineMode()
    el_tracker.sendCommand('clear_screen 0')
    
    # Print trial number on eyelink host monitor and output console
    status_msg = 'TRIAL number %d' % trial_dict['trial_num']
    el_tracker.sendMessage('TRIALID %d' % trial_dict['trial_num'])

    # Send status message to host PC
    el_tracker.sendCommand("record_status_message '%s'" % status_msg)

    # put tracker in idle/offline mode before recording
    el_tracker.setOfflineMode()
    
    # Start recording
    try:
        el_tracker.startRecording(1, 1, 1, 1) # arguments: sample_to_file, events_to_file, sample_over_link, event_over_link (1-yes, 0-no)
    except RuntimeError as error:
        print("ERROR:", error)
    
    # Allocate time for the tracker to cache some samples
    pylink.pumpDelay(100) 
    
    # Get eye used
    eye_used = get_eye_used(el_tracker)
    if eye_used is None and not EYETRACKER_OFF:
        print(f"Could not get eye used on trial {trial_dict['trial_num']}.")
        
    # Mark trial if eyetracker disconnected
    error = el_tracker.isRecording()
    if error is not pylink.TRIAL_OK:
        el_tracker.sendMessage('tracker_disconnected')
        print("Tracker disconnected.")
        
    # ------------------------------------------------------------------------
    # Start the trial loop
    while globalClock.getTime() < trial_end:
        update_target_onset = False
        record_rsvp_onset = False
        t = globalClock.getTime()
        
        # RSVP stream presentation
        if rsvp_idx < len(rsvp_sequence) and t >= next_pokemon_onset:
            record_rsvp_onset = True
            current_pokemon = pokemon_dict[rsvp_sequence[rsvp_idx]]
            rsvp_idx +=1
            next_pokemon_onset += RATE
            rsvpExp.nextEntry()
            # Gaze tracking until end of trial
            if not EYETRACKER_OFF:
                if not is_gaze_within_bounds(el_tracker, eye_used):
                    print(f"Trial {trial_dict['trial_num']}: Fixation broken")
                    thisExp.addData('fix.broken', 'fix.broken')
            
            if attention_cond == "FIX" and current_pokemon.name == target_pokemon:
                win.callOnFlip(lambda: rsvpExp.addData('target.onset', globalClock.getTime()))
                update_target_onset = True
                
        if current_pokemon is not None:
            current_pokemon.draw()
            
        # Peripheral stim logic
        if is_sim_trial:
            for pstim in pstim_to_draw:
                if trial_start + pstim_start_time <= t < trial_start + pstim_start_time + PERIPH_STIM_DURATION:
                    if attention_cond == "COV" and np.allclose(pstim.fillColor, PERIPHERAL_STIM_COLORS[TARGET_COLOR]):
                        update_target_onset = True
                    angle = pstim_info[pstim.name]['angle']
                    base_pos = pstim_info[pstim.name]['base_pos']
                    phase = pstim_info[pstim.name]['phase']
                    if angle is not None:
                        elapsed = t - pstim_start_time
                        angle_rad = np.deg2rad(angle)
                        offset = AMPLITUDE * np.sin(2*np.pi*FREQUENCY*elapsed + phase)
                        dy = offset * np.sin(angle_rad)
                        dx = offset * np.cos(angle_rad)
                        
                        pstim.pos = base_pos + (dx,dy)
                        
                    pstim.draw()
                    
        if is_seq_trial and pstim_idx < len(pstim_to_draw):
            current_pstim = pstim_to_draw[pstim_idx]
            onset_time = trial_start + pstim_onsets[pstim_idx]
            if onset_time <= t < onset_time + PERIPH_STIM_DURATION:
                if attention_cond == "COV" and np.allclose(pstim.fillColor, PERIPHERAL_STIM_COLORS[TARGET_COLOR]):
                    update_target_onset = True
                angle = pstim_info[current_pstim.name]['angle']
                base_pos = pstim_info[current_pstim.name]['base_pos']
                phase = pstim_info[current_pstim.name]['phase']
                if angle is not None:
                    elapsed = t - onset_time
                    angle_rad = np.deg2rad(angle)
                    offset = AMPLITUDE * np.sin(2*np.pi*FREQUENCY*elapsed + phase)
                    dx = offset * np.cos(angle_rad)
                    dy = offset * np.sin(angle_rad)
                    
                    current_pstim.pos = base_pos + (dx,dy)
                    
                current_pstim.draw()
                
            elif t >= onset_time + PERIPH_STIM_DURATION:
                pstim_idx +=1
                
        # record stimulus onsets and offsets when the window flips
        def on_flip(): 
            flip_time = globalClock.getTime(format='float')
            
            nonlocal last_target_onset, target_onset_recorded
            
            if is_sim_trial: 
                if trial_start + pstim_start_time <= t < trial_start + pstim_start_time + PERIPH_STIM_DURATION:
                    for pstim in pstim_to_draw:
                        if not pstim_onset_recorded_dict[pstim.name]:
                            thisExp.addData(f'{pstim.name}.onset', flip_time)
                            pstim_onset_recorded_dict[pstim.name] = True
                elif t >= trial_start + pstim_start_time + PERIPH_STIM_DURATION:
                    for pstim in pstim_to_draw:
                        if not pstim_offset_recorded_dict[pstim.name]:
                            thisExp.addData(f'{pstim.name}.offset', flip_time)
                            pstim_offset_recorded_dict[pstim.name] = True
                
            if is_seq_trial:
                if not pstim_onset_recorded_dict[current_pstim.name]:
                    thisExp.addData(f'{current_pstim.name}.onset', flip_time)
                    pstim_onset_recorded_dict[current_pstim.name] = True
                elif t >= onset_time + PERIPH_STIM_DURATION:
                    thisExp.addData(f'{current_pstim.name}.offset', flip_time)
                    pstim_offset_recorded_dict[current_pstim.name] = True
                    
            if record_rsvp_onset:
                rsvpExp.addData('run', run)
                rsvpExp.addData('block', trial_dict['presentation_cond'])
                rsvpExp.addData('trial', trial_dict['trial_num'])
                rsvpExp.addData('attention_cond', attention_cond)
                rsvpExp.addData('stim', current_pokemon.name)
                rsvpExp.addData('stim.onset', flip_time)
                
            rsvpExp.addData('stim.offset', flip_time)
                
            if update_target_onset and not target_onset_recorded:
                last_target_onset = flip_time
                target_onset_recorded = True
                
        win.callOnFlip(on_flip)
        win.flip()
        
        keys = kb.getKeys(keyList=[RESPONSE_KEY, 'escape'], waitRelease=False, clear=True)
        for key in keys:
            if key.name == 'escape':
                terminate_task()
            elif key.name == RESPONSE_KEY:
                press_time = key.rt
                rsvpExp.addData('press_time', press_time)
                press_times.append(press_time)
                if last_target_onset is not None:
                    rt = press_time - last_target_onset
                    rsvpExp.addData('rt', rt)
                    rts.append(rt)
                    if 0 < rt <= RESPONSE_WINDOW:
                        hit = 1
                        rsvpExp.addData('hit', hit)
    
    # Save offsets for pstims displayed in last 1 second of trial
    if is_seq_trial:
        thisExp.addData(f'{current_pstim.name}.offset', t)
    if is_sim_trial and pstim_start_time == 3:
        for pstim in pstim_to_draw:
            thisExp.addData(f'{pstim.name}.offset', t)
        
    # Save trial outputs to data file
    thisExp.addData('target_shown', target_shown)
    thisExp.addData('target.onset', last_target_onset) if target_shown else thisExp.addData('target.onset', '')
    thisExp.addData('press_times', press_times)
    thisExp.addData('rts', rts)
    thisExp.addData('hit', hit)
    thisExp.addData('keypresses', len(press_times)) # how many times pt responded in the trial
    thisExp.addData('trial.end', t)
        

    # Print trial data
    print(f"Trial block: {trial_dict['block_num']}, Trial: {trial_dict['trial_num']}, Target shown: {target_shown}, Keypresses: {len(press_times)}")

    thisExp.nextEntry()

    return last_target_onset

def practice_run(feat_cond, attention_cond, run, blanks_rsvps, trial_grids, trial_rsvps, sim_onsets):

    # Reset last target onset
    last_target_onset = None
    
    # Extract visuals for this run
    trial_dicts = create_trial_dicts(trial_grids, trial_rsvps, sim_onsets, True) # create trial dicts for this run
    
    # Do a drift check before the blank block
    drift_check()
    
    # Run blank block before trials
    thisExp.addData('prac_blank_block.start', globalClock.getTime(format='float'))
    last_target_onset = run_blank_block(blanks_rsvps[0], run, 1, attention_cond, feat_cond, last_target_onset, True)
    thisExp.addData('prac_blank_block.end', globalClock.getTime(format='float'))
    thisExp.nextEntry()

    # Run 1 SIM block and 1 SEQ block
    for trial_dict in trial_dicts:
        last_target_onset = run_trial(feat_cond, run, trial_dict, attention_cond, last_target_onset, True)
    rsvpExp.nextEntry()
    
    # Run blank block after trials
    thisExp.addData('prac_blank_block.start', globalClock.getTime(format='float'))
    last_target_onset = run_blank_block(blanks_rsvps[1], run, 2, attention_cond, feat_cond, last_target_onset, True)
    thisExp.addData('prac_blank_block.end', globalClock.getTime(format='float'))
    thisExp.nextEntry()
    
    # Feedback screen
    show_feedback(feat_cond, attention_cond)

def experiment_run(feat_cond, attention_cond, run, blanks_rsvps, trial_grids, trial_rsvps, sim_onsets):

    # Reset last target onset
    last_target_onset = None
    
    # Extract visuals for this run
    trial_dicts = create_trial_dicts(trial_grids, trial_rsvps, sim_onsets) # create trial dicts for this run
    
    # Do a drift check before the blank block
    drift_check()

    # Run blank block before trials
    thisExp.addData('blank_block.start', globalClock.getTime(format='float'))
    last_target_onset = run_blank_block(blanks_rsvps[0], run, 1, attention_cond, feat_cond, last_target_onset)
    thisExp.addData('blank_block.end', globalClock.getTime(format='float'))
    thisExp.nextEntry()

    # Run all trials using trial dictionaries, updating accuracy
    for trial_dict in trial_dicts:
        last_target_onset = run_trial(feat_cond, run, trial_dict, attention_cond, last_target_onset)
    rsvpExp.nextEntry()
    
    # Run blank block after trials
    thisExp.addData('blank_block.start', globalClock.getTime(format='float'))
    last_target_onset = run_blank_block(blanks_rsvps[1], run, 2, attention_cond, feat_cond, last_target_onset)
    thisExp.addData('blank_block.end', globalClock.getTime(format='float'))
    thisExp.nextEntry()
    
    # Feedback screen
    show_feedback(feat_cond, attention_cond)
  
###### WELCOME SCREEN ##################################################################################################################

# Clear the window and set start time
win.flip() 
thisExp.addData('welcome.start', globalClock.getTime(format='float'))
thisExp.status = STARTED

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

# Calibrate eyetracker if on
if not EYETRACKER_OFF:
    thisExp.addData('calibration.start', globalClock.getTime(format='float'))
    try:
        el_tracker.doTrackerSetup()
    except RuntimeError as err:
        print('ERROR:', err)
        el_tracker.exitCalibration()
    thisExp.addData('calibration.end', globalClock.getTime(format='float'))
else: # Wait for space or escape key
    keys = event.waitKeys(keyList=['space', 'escape'])
    if 'escape' in keys:
        terminate_task()

###### EXPERIMENT #####################################################################################################################

# Clear the window
win.flip() 

# Generate or load the RSVP sequences and peripheral stimulus grids for all runs
prac_blank_rsvps, exp_blank_rsvps = generate_blank_rsvps()
prac_trial_rsvps, exp_trial_rsvps = generate_trial_rsvps()
prac_sim_onsets, exp_sim_onsets = generate_sim_onsets()
prac_trial_grids, exp_trial_grids = assign_grids(feat_cond)

# FIX 
if any(r in run_list for r in [1, 2, 3]):
    attention_cond = 'FIX'
    show_instructions(feat_cond, attention_cond)
    for run in range(1, PRACTICE_RUNS+1):
        practice_run(feat_cond, attention_cond, run, prac_blank_rsvps[run], prac_trial_grids[run], prac_trial_rsvps[run], prac_sim_onsets[run])
    ready_screen()
    for run in run_list:
        if run <= RUNS_PER_COND:
            experiment_run(feat_cond, attention_cond, run, exp_blank_rsvps[run], exp_trial_grids[run], exp_trial_rsvps[run], exp_sim_onsets[run])

# COV
if any(r in run_list for r in [4, 5, 6]):
    attention_cond = 'COV'
    show_instructions(feat_cond, attention_cond)
    for run in range(1, PRACTICE_RUNS+1):
        practice_run(feat_cond, attention_cond, run, prac_blank_rsvps[run], prac_trial_grids[run], prac_trial_rsvps[run], prac_sim_onsets[run])
    ready_screen()
    for run in run_list:
        if run >= RUNS_PER_COND:
            experiment_run(feat_cond, attention_cond, run, exp_blank_rsvps[run], exp_trial_grids[run], exp_trial_rsvps[run], exp_sim_onsets[run])

###### END EXPERIMENT ##################################################################################################################

# Draw thank you text 
thanks_text = visual.TextStim(win=win, text="Thanks for coming to the Pokémon Party!", font='Arial', units='deg', pos=(0, 0), height=1.2, wrapWidth=1700, 
    color='black', colorSpace=CLR_SPC)
thanks_text.draw()
win.flip()

# Close experiment window and save data when space is pressed
keys = event.waitKeys(keyList=['space'])
terminate_task() 