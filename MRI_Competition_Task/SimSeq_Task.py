from psychopy import prefs, plugins, sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware, monitors
from psychopy.constants import (NOT_STARTED, STARTED, STOPPED, FINISHED, PRESSED, RELEASED, FOREVER, priority)
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy.hardware import keyboard
import math
import numpy as np
import itertools
import random
import os
import sys
import time
logging.console.setLevel(logging.ERROR)

###### EXPERIMENT PARAMETERS #######################################################################################################

# Initialize the global clock
globalClock = core.Clock()

# Dictionaries of feature conditions
"""Color condition: 2 attention conditions, 2 presentation conditions
                    6 runs total, 2 blank blocks per run, 12 blocks of sim/seq per run, 3 trials per block, 
                    predetermined block design"""
color_condition = {'ATTENTION_CONDS': ['FIX', 'COV'], 
                   'PRESENTATION_CONDS': ['SIM', 'SEQ'],
                   'NUM_RUNS': 6, 
                   'NUM_BLANKS': 2,
                   'NUM_BLOCKS': 12, 
                   'NUM_TRIALS': 3,
                   'BLOCK_DESIGN': [('RVF','SIM'),('LVF','SEQ'),('RVF','SEQ'),('LVF','SIM'),('RVF','SEQ'),('LVF','SIM'),
                                    ('RVF','SIM'),('LVF','SEQ'),('RVF','SIM'),('LVF','SEQ'),('RVF','SEQ'),('LVF','SIM')]}

# Size
PERIPHERAL_STIM_SIZE = 1.75 #DVA; size of each peripheral stimulus (circle)
POKEMON_SIZE = [1.5, 1.5] # DVA, size of the pokemon during RSVP
NUM_PSTIMS = 4 # number of peripheral stimuli
GRID_SIZE = 4 #DVA; height and width of the peripheral stimulus grid
ECCENTRICITY = 7 #DVA from the center of the grid to the center of each peripheral stimulus
PERIPHERAL_STIM_COLORS = ['red', 'blue', 'green', 'yellow', 'magenta', 'cyan'] 

# Timing
BLANK_BLOCK_DURATION = 16 # seconds
PERIPH_STIM_DURATION = 1 # seconds
ISI = 0.033 # seconds; blank fixation between peripheral stimuli for SEQ condition
TRIAL_DURATION = PERIPH_STIM_DURATION*NUM_PSTIMS # seconds, not including ISIs
PSTIM_TARGET_FREQ = [1,3] # pstim color targets will occur every 1-3 trials
POKEMON_TARGET_FREQ = [3.75,7.5] # pokemon targets will occur every 3.75-7.5s in the RSVP

# Response Parameters
RSVP_RATE = 0.25 # seconds; new pokemon every 250ms in the RSVP
RESPONSE_KEY = 'space'

###### EXPERIMENT SETUP #######################################################################################################

exp_name = 'SimSeq-DV'
exp_info = {'Participant ID': '9999', 
            'Session': '001',
            }
dlg = gui.DlgFromDict(dictionary=exp_info, title=exp_name)
if dlg.OK == False:
    core.quit()  

# Establish data output directory
time_str = time.strftime("_%m_%d_%Y", time.localtime())
root_dir = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(root_dir, 'data', f"{exp_name}_{exp_info['Participant ID']}_Session{exp_info['Session']}_{time_str}")
os.makedirs(data_folder, exist_ok=True)
filename = os.path.join(data_folder, f"{exp_name}_{exp_info['Participant ID']}_Session{exp_info['Session']}")

# Window setup (will need to be adjusted to match the MRI monitor)
win = visual.Window(fullscr=True,color=[0,0,0], screen=0, 
                    size = [1512,982], monitor='testMonitor',
                    winType='pyglet', allowStencil=False,
                    blendMode='avg', useFBO=False,
                    colorSpace='rgb', units='deg')

# Create an experiment handler to manage the data file
thisExp = data.ExperimentHandler(name=exp_name, version='', extraInfo=exp_info,
                                runtimeInfo=None, originPath=os.path.abspath(__file__),
                                savePickle=True, saveWideText=True,
                                dataFileName=filename)

###### INITIALIZE VISUAL COMPONENTS #######################################################################################################

kb = keyboard.Keyboard()

# Load Pokémon images from folder
pokemon_dir = os.path.join(root_dir, 'pokemon_lightgray')
pokemon_names = ["Bulbasaur", "Pikachu", "Squirtle", "Charmander", "Magikarp", "Raticate", "Pidgey",
    "Metapod", "Jigglypuff", "Butterfree", "Psyduck", "Caterpie", "Krabby",
    "Haunter", "Vulpix", "Eevee", "Sandshrew", "Wartortle", "Charmeleon", "Clefairy",
    "Ponyta", "Mankey"]
pokemon_dict = {name: visual.ImageStim(win, name=name, image=os.path.join(pokemon_dir, f"{i+1:03}.png"))
    for i, name in enumerate(pokemon_names)} # dictionary of the pokemon where the key is their name and the values are the ImageStim

# Components for welcome screen
welcome_text = visual.TextStim(win, pos=(0,0), height= 1.5, wrapWidth=27, text=("Welcome to the Pokémon Party game!"))
# Subset of pokemon to be shown on the welcome screen with their position and size
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

# Components for instructions screen
instructions_text = visual.TextStim(win, pos = (0,0), wrapWidth=27, text=())

# Components for practice and experiment blocks
pstim_dict = {color: visual.Circle(win, name=color, pos = (), size = 2, edges = 32, color = color, units = 'deg', anchor = 'center')
          for color in PERIPHERAL_STIM_COLORS} # dictionary of the pstim where the key is their color and the values are the Circle

end_text = visual.TextStim(win, wrapWidth=27, text=())

# Calculate the coordinates of the center of the grid based on the given parameters
cent2cent_spacing = GRID_SIZE - PERIPHERAL_STIM_SIZE # 2.25DVA; distance from center to center of peripheral stimuli
offset = cent2cent_spacing / 2 # 1.125DVA; how much to move in x and y from the center of the grid to the center of each peripheral stimulus
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

# Components for the end of the experiment
thanks_text = visual.TextStim(win, wrapWidth=30, text=("Thanks for coming to the Pokémon Party!"))

###### FUNCTIONS #######################################################################################################

def end_task():
    """ Saves data and closes the window."""
    thisExp.nextEntry()
    thisExp.addData('experiment.stopped', globalClock.getTime(format='float'))
    
    print("Task ended.")
    thisExp.saveAsWideText(filename + '.csv', delim='auto')
    thisExp.saveAsPickle(filename)
    logging.flush()
    
    if win is not None:
        win.clearAutoDraw()
        win.flip()

    thisExp.abort() # or data files will save again on exit
    win.close()
    core.quit()
    sys.exit()

def draw_comp(comp, t, tThisFlipGlobal, frameN):
    """ Draws the stimuli during the trials. """
    
    comp.tStart = t
    comp.tStartRefresh = tThisFlipGlobal
    comp.frameNStart = frameN
    comp.status = STARTED
    if not isinstance(comp, keyboard.Keyboard):
        thisExp.addData(f'{comp.name}.started', comp.tStartRefresh)
        comp.setAutoDraw(True)
    else:
        thisExp.addData('kb.started', comp.tStartRefresh)

def erase_comp(comp, t, tThisFlipGlobal, frameN):
    """ Erases the stimuli displayed in trials. """
    comp.tStop = t
    comp.tStopRefresh = tThisFlipGlobal
    comp.frameNStop = frameN
    comp.status = FINISHED
    if not isinstance(comp, keyboard.Keyboard):
        comp.setAutoDraw(False)
    else:
        thisExp.addData('kb.stopped', comp.tStopRefresh)

def generate_all_visuals(feature_condition, target_pokemon, target_color):
    """ Generates all visual stimuli for running the feature condition (i.e., the RSVP Pokémon sequences and peripheral stimulus grids).
    
    NOTES:
        -Same target pokémon and target peripheral stim color used for all runs in the feature condition.
        -Target pokémon presentation is based on time interval between targets (s); target pstim presentation is based on number of trials between targets.
        -Each trial in a run will aim to have a unique peripheral stimulus grid (unless all are taken, then resampled).
        -The visuals for the run are repeated across attention conditions (e.g., run 1 and run 4 of the color condition will look exactly the same). 
        
    Returns: 
            all_run_visuals (dict): A dictionary where the keys are the runs and values are dictionaries of parameters for each trial in the run. 
                Each trial dictionary will have the following information:
                    'block_num': 1,
                    'trial_num': 1,
                    'vf': 'right',
                    'presentation_cond': 'SIM',
                    'periph_stim_grid': ('red', 'blue', 'green', 'yellow'),
                    'rsvp_sequence': ['Bulbasaur', 'Squirtle', 'Pikachu', ..., 'Meowth']
            blanks_rsvps (dict): A dictionary where the keys are the runs and the values are 2 lists of RSVPs composed of distractor pokémon.
    """
    num_blocks = feature_condition['NUM_BLOCKS']
    num_trials = feature_condition['NUM_TRIALS']
    block_design = feature_condition['BLOCK_DESIGN']
    num_runs = feature_condition['NUM_RUNS']
    num_blanks = feature_condition['NUM_BLANKS']
    attention_conds = feature_condition['ATTENTION_CONDS']
    
    num_attention_conds = len(attention_conds)
    num_unique_runs = num_runs//num_attention_conds # visuals are repeated across the attention conds 
    trials_per_run = num_blocks*num_trials
    total_run_duration = trials_per_run*TRIAL_DURATION # does not include ISIs or blank blocks

    # ---------------- Generate blank block RSVPs for all runs ----------------

    num_blank_rsvps = int(num_blanks*num_runs//2) # total number of unique RSVP sequences for the blank blocks (each repeated for each attention_cond)
    num_pokemon_needed = int(BLANK_BLOCK_DURATION//RSVP_RATE) # how many pokemon needed for one blank block RSVP
    
    blank_sequences = []
    for block in range(num_blank_rsvps):
        target_indices = [] # list of time t when target pokemon appears in the blank block, with t=0 being the start of the block
        t_pokemon = random.uniform(*POKEMON_TARGET_FREQ) # randomly pick the first target pokemon onset
        while t_pokemon < BLANK_BLOCK_DURATION:
            target_idx = int(t_pokemon//RSVP_RATE)
            if target_idx < num_pokemon_needed:
                target_indices.append(target_idx)
            t_pokemon += random.uniform(*POKEMON_TARGET_FREQ) # add a random time interval between last target and next target
            
        rsvp_sequence = []
        for idx in range(num_pokemon_needed):
            if idx in target_indices:
                rsvp_sequence.append(target_pokemon)
            else:
                distractor = random.choice([p for p in pokemon_names if p != target_pokemon])
                # ensure that there are no back to back repetitions of pokemon in the sequence
                if idx !=0 and distractor == rsvp_sequence[-1]:
                    distractor = random.choice([p for p in pokemon_names if p != target_pokemon and p != rsvp_sequence[-1]])
                rsvp_sequence.append(distractor)
        blank_sequences.append(rsvp_sequence) # blank_sequences is a list of 6 inner lists of pokemon sequences

    blank_sequences *= num_attention_conds # repeat the unique sequences for each attention condition

    # in a dictionary, assign one rsvp to start blank block and one end blank block to each run in the feature condition
    blanks_rsvps = {}
    seq_idx = 0
    for run in range(1, num_runs+1):
        blanks_rsvps[run] = [blank_sequences[seq_idx], blank_sequences[seq_idx+1]]
        seq_idx += num_attention_conds
        
    # ---------------- Peripheral Stimulus Grid Setup ----------------
    
    # Create list of all possible 4-color permutations
    color_combos = itertools.combinations(PERIPHERAL_STIM_COLORS, NUM_PSTIMS) # all the different ways to create a subset of NUM_PSTIMS colors from 6
    all_grids = []
    for combo in color_combos:
        permutations = itertools.permutations(combo) # all the different ways to order those subsets of NUM_PSTIMS colors 
        all_grids.extend(permutations)
    random.shuffle(all_grids)

    # Split grids into those with and without the target color
    grids_with_target = [grid for grid in all_grids if target_color in grid]
    grids_without_target = [grid for grid in all_grids if target_color not in grid]

    # ---------------- Function to create trial dicts for a specific run ----------------
    def create_trial_dicts(run_number):
        
        # Trials that will have target pstim
        target_pstim_trials = [] # list of trial indices 
        pstim_trial_idx = random.randint(0,PSTIM_TARGET_FREQ[1]) # randomly pick the first target trial 
        while pstim_trial_idx < trials_per_run:
            target_pstim_trials.append(pstim_trial_idx)
            pstim_trial_idx += random.randint(*PSTIM_TARGET_FREQ) # randomly choose how many trials before the next target color onset
            
        # Assign grids to trials
        grid_assignments = [None] * trials_per_run
        used_grids = set()
        for trial_idx in range(trials_per_run):
            if trial_idx in target_pstim_trials:
                available = [grid for grid in grids_with_target if grid not in used_grids] or grids_with_target # allow reuse if all lists with target have been used
            else:
                available = [grid for grid in grids_without_target if grid not in used_grids] or grid_without_target
            grid = random.choice(available)
            grid_assignments[trial_idx] = grid # list of all the grids to be used for the entire run
            used_grids.add(grid) # remove grid from pool of available grids, each trial will have a unique grid 
        
        # ---------------- Pokemon RSVP List Setup for Experiment Trials ----------------
        
        # Determine which trials will have the target pokemon
        target_pokemon_times = [] # list of time t when target pokemon appears in the run, with t=0 being the start of the run
        t_pokemon = random.uniform(0, POKEMON_TARGET_FREQ[1]) # randomly pick the first target pokemon onset 
        while t_pokemon < total_run_duration:
            target_pokemon_times.append(t_pokemon)
            t_pokemon += random.uniform(*POKEMON_TARGET_FREQ) # add a random time interval between last target and next target

        # Map pokemon target onsets to trials
        target_onsets = {i: [] for i in range(trials_per_run)} # dictionary of target pokemon onsets within the trial, with t=0 being start of the trial
        for t in target_pokemon_times:
            trial_idx = int(t // TRIAL_DURATION)
            if trial_idx < trials_per_run:
                onset_within_trial = t - (trial_idx * TRIAL_DURATION)
                target_onsets[trial_idx].append(onset_within_trial) # add target pokemon onset to the dictionary, the keys are trial indices and the values are the onsets

        # ---------------- Build Trial Dictionaries ----------------
        trial_dicts = []
        for block_idx, (vf, presentation_cond) in enumerate(block_design):
            for trial in range(num_trials):
                trial_idx = block_idx * num_trials + trial
                assigned_grid = grid_assignments[trial_idx]
                
                # Create RSVP sequence
                total_trial_pokemon = int(TRIAL_DURATION // RSVP_RATE) # total number of pokemon seen each trial (not including ISI)
                rsvp_sequence = []
                for i in range(total_trial_pokemon):
                    timepoint = i * RSVP_RATE
                    if any(abs(timepoint - t) < RSVP_RATE / 2 for t in target_onsets[trial_idx]):
                        rsvp_sequence.append(target_pokemon)
                    else:
                        distractor = random.choice([p for p in pokemon_names if p != target_pokemon])
                        # ensure that there are no back to back repetitions of pokemon in the sequence
                        if i !=0 and distractor == rsvp_sequence[-1]:
                            distractor = random.choice([p for p in pokemon_names if p != target_pokemon and p != rsvp_sequence[-1]])
                        rsvp_sequence.append(distractor)
                
                trial_dicts.append({
                    'block_num': block_idx + 1,
                    'trial_num': trial + 1,
                    'vf': vf,
                    'presentation_cond': presentation_cond,
                    'periph_stim_grid': assigned_grid,
                    'rsvp_sequence': rsvp_sequence
                })
                
        return trial_dicts

    # ---------------- Create trial dicts for all the runs in the feature condition ----------------
    all_run_visuals = {}
    unique_run_visuals = [create_trial_dicts(i+1) for i in range(num_unique_runs)]
    for run_num in range(1, num_runs+1):
        all_run_visuals[f"run{run_num}"] = unique_run_visuals[(run_num-1)%num_unique_runs]

    return all_run_visuals, blanks_rsvps

def run_trial(trial_dict, attention_cond, target_pokemon, target_color):
    """
    Function to run a single trial.

    Args:
       trial_dict(dict): Contains 'block_num', 'trial_num', 'vf', 'presentation_cond', 'periph_stim_grid','rsvp_sequence'
       attention_cond (str): 'FIX' (looking for a pokemon) or 'COV' (looking for a colored circle)
       target_pokemon (str): Name of pokemon to look for
       target_color (str): Name of color to look for
    """

    # Reset clock and frame
    routineTimer.reset()
    frameN = -1

    # Reset indices
    rsvp_index = 0
    pstim_index = 0
    pstim_onset = 0 # for SEQ trials
    target_pstim_onset = None
    target_pokemon_onset = None

    # Save block and trial number to the data file
    thisExp.addData('pres_cond', trial_dict['presentation_cond'])
    thisExp.addData('block', trial_dict['block_num'] )
    thisExp.addData('trial', trial_dict['trial_num'])

     # Get trial parameters from the trial dictionary
    vf = trial_dict['vf']
    is_seq_trial = trial_dict['presentation_cond'] == 'SEQ'
    is_sim_trial = trial_dict['presentation_cond'] == 'SIM'
    pstim_grid = trial_dict['periph_stim_grid']
    rsvp_sequence = trial_dict['rsvp_sequence']

    # Save peripheral stim grid to data file
    thisExp.addData('periph_stim', pstim_grid)

    # Reset lists of key presses, rt, and peripheral stim to draw. Set all trial visual component status to not started
    continueRoutine = True
    kb.clearEvents()
    kb.rt = []
    kb.keys = []
    kb_allKeys = []
    pstim_to_draw = []

    components = [pstim_dict, pokemon_dict, kb]
    for comp in components:
        if isinstance(comp, dict):
            for i in comp:
                comp[i].tStart = None
                comp[i].tStop = None
                comp[i].tStopRefresh = None
                if hasattr(comp[i], 'status'):
                    comp[i].status = NOT_STARTED
        else:
            comp.tStart = None
            comp.tStop = None
            comp.tStartRefresh = None
            comp.tStopRefresh = None
            if hasattr(comp, 'status'):
                comp.status = NOT_STARTED

    # Extract peripheral stim colors and locations
    if vf == 'RVF':
        grid_positions = [rvf_topright, rvf_topleft, rvf_botleft, rvf_botright]
    elif vf == 'LVF':
        grid_positions = [lvf_topright, lvf_topleft, lvf_botleft, lvf_botright]

    # Map colors to their positions based on their order in pstim_grid
    for color, pos in zip(pstim_grid, grid_positions):
        pstim = pstim_dict[color]
        pstim.pos = pos
        pstim_to_draw.append(pstim) # add all the pstims to be drawn to a list

    # Used by SIM only: select random start time for peripheral stimuli
    pstim_start_time = random.choice([0, 1, 2, 3])

    # For SEQ only: Determine how many more pokemon need to be added to the end of the RSVP sequence to account for additional ISI time
    if is_seq_trial:
        additional_pokemon = math.ceil(NUM_PSTIMS * ISI / RSVP_RATE) 
        distractors = [pokemon_name for pokemon_name in pokemon_names if pokemon_name != target_pokemon] # only non-target pokemon will be added to the RSVP
        rsvp_sequence.extend(random.choices(distractors, k=additional_pokemon))

    # Save RSVP sequence to the data file
    thisExp.addData('rsvp_stim', rsvp_sequence)

    # Start the trial loop
    while continueRoutine:
        t = routineTimer.getTime() # t = 0 is start of the trial
        tThisFlip = win.getFutureFlipTime(clock=routineTimer) # when the next screen flip will be, relative to the start of the trial
        tThisFlipGlobal = win.getFutureFlipTime(clock=globalClock) # when the next screen flip will be, relative to the start of the experiment
        frameN += 1

        if rsvp_index < len(rsvp_sequence): # make sure we don't go past the last pokemon in the rsvp
            current_pokemon = pokemon_dict[rsvp_sequence[rsvp_index]]
            is_final_pokemon = rsvp_index == len(rsvp_sequence) - 1 # check if this is the last pokemon in the RSVP

            # Visual flow for SIM condition
            if is_sim_trial:
                # Draw pokemon RSVP stream
                if current_pokemon.status == NOT_STARTED and tThisFlip >= 0: # Start RSVP stream at the start of the trial
                    draw_comp(current_pokemon, t, tThisFlipGlobal, frameN)
                    current_pokemon_onset = current_pokemon.tStartRefresh # Record time when the pokemon was first presented
                    if rsvp_sequence[rsvp_index] == target_pokemon and target_pokemon_onset == None:
                        target_pokemon_onset = current_pokemon_onset # If pokemon is target pokemon, save first onset time

                if current_pokemon.status == STARTED and tThisFlipGlobal > current_pokemon_onset + RSVP_RATE: # erase the current pokemon after rsvp rate elapses
                    erase_comp(current_pokemon, t, tThisFlipGlobal, frameN)
                    if not is_final_pokemon:
                        rsvp_index += 1 # move to the next pokemon in the sequence unless this is the last one
                        pokemon_dict[rsvp_sequence[rsvp_index]].status = NOT_STARTED # reset status for the next pokemon to be drawn in case it has already been shown
                
                # Draw the peripheral stimulus grid all at once at a random timepoint
                for pstim in pstim_to_draw:
                    if pstim.status == NOT_STARTED and tThisFlip >= pstim_start_time:
                        draw_comp(pstim, t, tThisFlipGlobal, frameN)
                        if pstim.color == target_color and target_pstim_onset == None:
                            target_pstim_onset = pstim.tStartRefresh # If target color in grid, save onset time

                    if pstim.status == STARTED and tThisFlipGlobal > pstim.tStartRefresh + PERIPH_STIM_DURATION: # Erase grid after pstim duration
                        erase_comp(pstim, t, tThisFlipGlobal, frameN)

            # Visual flow for SEQ condition
            if is_seq_trial:
                # Draw pokemon RSVP stream
                pokemon_duration = 0.132 if is_final_pokemon else RSVP_RATE # last pokemon is displayed for 132ms in SEQ to account for ISI

                if current_pokemon.status == NOT_STARTED and tThisFlip >= 0: # Start RSVP stream at the start of the trial
                    draw_comp(current_pokemon, t, tThisFlipGlobal, frameN)
                    current_pokemon_onset = current_pokemon.tStartRefresh # Record time when the pokemon was first presented
                    if rsvp_sequence[rsvp_index] == target_pokemon and target_pokemon_onset == None:
                        target_pokemon_onset = current_pokemon_onset # If pokemon is target pokemon, save first onset time

                if current_pokemon.status == STARTED and tThisFlipGlobal > current_pokemon_onset + pokemon_duration: # erase the current pokemon after rsvp rate elapses
                    erase_comp(current_pokemon, t, tThisFlipGlobal, frameN)
                    if not is_final_pokemon:
                        rsvp_index += 1 # move to the next pokemon in the sequence unless this is the last one
                        pokemon_dict[rsvp_sequence[rsvp_index]].status = NOT_STARTED # reset status for the next pokemon to be drawn

                # Draw peripheral stimuli one at a time
                if pstim_index < len(pstim_to_draw):
                    current_pstim = pstim_to_draw[pstim_index]

                    if current_pstim.status == NOT_STARTED and tThisFlip >= pstim_onset:
                        draw_comp(current_pstim, t, tThisFlipGlobal, frameN)
                        if current_pstim.color == target_color and target_pstim_onset == None:
                            target_pstim_onset = current_pstim.tStartRefresh # If target color shown, save onset time

                    if current_pstim.status == STARTED and tThisFlipGlobal > current_pstim.tStartRefresh + PERIPH_STIM_DURATION:
                        erase_comp(current_pstim, t, tThisFlipGlobal, frameN)
                        pstim_index += 1
                        pstim_onset = current_pstim.tStopRefresh + ISI

        # Start checking for key presses at start of trial
        if kb.status == NOT_STARTED and tThisFlip >= 0:
            draw_comp(kb, t, tThisFlipGlobal, frameN)

        # Check keypresses once per frame
        elif kb.status == STARTED:
            keys = kb.getKeys(keyList=[RESPONSE_KEY, 'escape'], waitRelease=False)
            for key in keys:
                kb_allKeys.append(key)
                if key.name == 'escape':
                    end_task()
                if key.name == RESPONSE_KEY:
                    print("in the if loop")
                    if not kb.keys:
                        thisExp.timestampOnFlip(win, 'kb.pressed')
                        kb.keys = key.name
                        kb.rt = kb_allKeys[0].rt # get the reaction time of the first response keypress

        # End the trial after the total duration
        if (is_sim_trial and tThisFlip >= TRIAL_DURATION) or (is_seq_trial and tThisFlip >= TRIAL_DURATION + (NUM_PSTIMS * ISI)):
            continueRoutine = False 
            for comp in components: # make sure all components are erased
                if isinstance(comp, dict):
                    for i in comp:
                        erase_comp(comp[i], t, tThisFlipGlobal, frameN)
                else:
                    erase_comp(comp, t, tThisFlipGlobal, frameN)

        if continueRoutine:
            win.flip()  # refresh the screen

    # Initialize trial outcome
    trial_outcome = "undefined"
    accuracy = None  # default to None

    if attention_cond == 'FIX':
        target_shown = target_pokemon in rsvp_sequence # T or F
        target_onset = target_pokemon_onset

    elif attention_cond == 'COV':
        target_shown = target_color in pstim_grid # T or F
        target_onset = target_pstim_onset

    response_keys = [key for key in kb_allKeys if key.name == RESPONSE_KEY]
    first_press = response_keys[0] if response_keys else None

    # Calculate accuracy
    if target_shown: # case 1: target shown in trial
        if not response_keys: # no key pressed
            accuracy = 0
            trial_outcome = "miss"
        elif first_press and first_press.rt >= target_onset: # key pressed after target was displayed; ANDREW: max duration after target?
            print("target shown", "RT:", first_press.rt, "target onset time:", target_onset)
            accuracy = 1
            trial_outcome = "hit"
        elif first_press and first_press.rt < target_onset: # key pressed before target displayed
            print("target shown", "RT:", first_press.rt, "target onset time:", target_onset)
            accuracy = 0
            trial_outcome = "early press"
    else: # case 2: target not shown in trial
        if not response_keys: # no key pressed and target not displayed
            accuracy = 1
            trial_outcome = "correct rejection"
        else:
            accuracy = 0 # key pressed when no target displayed
            trial_outcome = "false alarm"

    # Save accuracy and reaction time to file
    thisExp.addData('target_shown', target_shown)
    thisExp.addData('accuracy', accuracy)
    thisExp.addData('rt', kb.rt)
    thisExp.addData('outcome', trial_outcome)

    # Print trial data
    print(f"Block {trial_dict['block_num']}, Trial {trial_dict['trial_num']}, Accuracy {accuracy}: {trial_outcome}")

    thisExp.nextEntry()

    return accuracy

def run_feature_cond(feature_condition, target_pokemon, target_color):
    """
    Runs all of the blocks in a single feature condition. 
    
    Args:
        feature_condition (dict): Dictionary of all parameters for the feature condition defined at start of script.
        target_pokemon (str): Name of target pokémon.
        target_color (str): Name of color.
    """

    runs_per_attncond = feature_condition['NUM_RUNS']//len(feature_condition['ATTENTION_CONDS'])
    
    # Step 1: Generate the RSVP sequences and peripheral stimulus grids for all runs in the feature condition
    all_run_visuals, blanks_rsvps = generate_all_visuals(feature_condition, target_pokemon, target_color)

    # Step 2: Perform each run
    for run in range(1, feature_condition['NUM_RUNS'] + 1):
        run_accuracy = 0
        run_idx = f"run{run}"
        run_visuals = all_run_visuals[run_idx] # trial dicts for all the trials in this run
        blank_block_rsvp = blanks_rsvps[run] # 2 RSVP lists for the blank blocks in this run
        
        # Show instructions if this is the first run in the attention condition; kind of hacky only works with 2 attn conds
        if run <= runs_per_attncond:
            attention_cond = feature_condition['ATTENTION_CONDS'][0]
            if run == 1:
                instructions_text.text = (
                "There's a Pokémon Party happening right now, and the Pokémon are playing hide and seek!\n\n"
                f"The Pokémon are having trouble finding {target_pokemon}! Can you help them?\n\n"
                f"Press the button as fast as you can every time you see {target_pokemon}.\n\n\n"
                "Ready to start playing?")
            
        elif run > runs_per_attncond:
            attention_cond = feature_condition['ATTENTION_CONDS'][1]
            if run == 1+runs_per_attncond:
                instructions_text.text = (
                "There's a Pokémon Party happening right now, and the Pokémon are getting hungry!\n\n"
                f"The Pokémon like to eat {target_color} circles! Can you help feed them?\n\n"
                f"Press the button as fast as you can every time you see a {target_color} circle.\n\n\n"
                "Ready to start playing?")
        instructions_text.draw()
        win.flip()
        
        # Save attention condition and instructions start time to data file
        thisExp.addData('attention_cond', attention_cond)
        thisExp.addData('instructions.started', globalClock.getTime(format='float'))
        
        # Wait for space or esc 
        keys = event.waitKeys(keyList=['space', 'escape'])
        if 'escape' in keys:
            end_task()
        
        # Save the instructions stop time to the data file and move to the next row
        thisExp.addData('instructions.stopped', globalClock.getTime(format='float')) 
        
        # Blank block; display an RSVP of distractors
        for pokemon in blank_block_rsvp[0]:
            pokemon_dict[pokemon].pos = (0,0)
            pokemon_dict[pokemon].draw()
            win.flip()
            core.wait(RSVP_RATE)
            keys = kb.getKeys(keyList=['escape'], waitRelease=False)
            if 'escape' in keys:
                end_task()
        
        thisExp.addData('blank_block.stopped', globalClock.getTime(format='float'))
        
        # Feed trial dictionary and attention condition to the run trial function
        thisExp.nextEntry() 
        for trial in run_visuals:
            run_accuracy += run_trial(trial, attention_cond, target_pokemon, target_color)
        
        # Blank block; display an RSVP of distractors 
        for pokemon in blank_block_rsvp[1]:
            pokemon_dict[pokemon].pos = (0,0)
            pokemon_dict[pokemon].draw()
            win.flip()
            core.wait(RSVP_RATE)
            keys = kb.getKeys(keyList=['escape'], waitRelease=False)
            if 'escape' in keys:
                end_task()
            
        thisExp.addData('blank_block.stopped', globalClock.getTime(format='float'))
        
        # Step 3: Calculate hit rate and show feedback statement at the end of the run based on attention condition
        hit_rate = run_accuracy / (feature_condition['NUM_RUNS']*feature_condition['NUM_TRIALS']) ## ANDREW: different feedback after the run depending on hit rate?
        if attention_cond == 'FIX':
            end_text.text = (f"Great job finding {target_pokemon}!")
            pokemon_dict[target_pokemon].pos = (0, -5) 
            pokemon_dict[target_pokemon].size = (5,5)
            pokemon_dict[target_pokemon].draw()
        elif attention_cond == 'COV':
            end_text.text = (f"Great job feeding the Pokémon {target_color} circles!")
            pstim_dict[target_color].pos = (0,-5)
            pokemon_dict[target_pokemon].size = (5,5)
            pstim_dict[target_color].draw()
        end_text.draw()
        win.flip()
        
        # Reset size and position of the target pokemon
        pokemon_dict[target_pokemon].size = POKEMON_SIZE
        pokemon_dict[target_pokemon].pos = (0,0)
        
        keys = event.waitKeys(keyList=['space', 'escape'])
        if 'escape' in keys:
            end_task()

###### WELCOME SCREEN #######################################################################################################

# Start experiment clock and clear the window
routineTimer = core.Clock()
win.flip() 

# Mark the experiment as started
exp_info['expStart'] = data.getDateStr(format = '%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6)
thisExp.status = STARTED

# Draw welcome screen with Pokémon images
win.mouseVisible = False  # Hide mouse cursor
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

# Wait for space or escape key
keys = event.waitKeys(keyList=['space', 'escape'])
if 'escape' in keys:
    end_task()

###### PRACTICE BLOCK #######################################################################################################

###### EXPERIMENT BLOCK #######################################################################################################

# Run each feature condition
run_feature_cond(color_condition, 'Pikachu', 'red')

###### END EXPERIMENT #######################################################################################################

# Draw thank you text and run function to save data
thanks_text.draw()
win.flip()
keys = event.waitKeys(keyList=['space'])

end_task() 