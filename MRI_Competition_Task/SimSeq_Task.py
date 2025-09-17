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
SIM_TRIAL_DURATION = PERIPH_STIM_DURATION*NUM_PSTIMS # seconds, not including ISIs
SEQ_TRIAL_DURATION = SIM_TRIAL_DURATION + (NUM_PSTIMS-1)*ISI # assumes no ISI between trials
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

def draw_comp(comp, t, tThisFlip, tThisFlipGlobal, frameN):
    """ Draws the stimuli during the trials. """
    comp.tStart = tThisFlip
    comp.tStartRefresh = tThisFlipGlobal
    comp.frameNStart = frameN
    comp.status = STARTED
    if not isinstance(comp, keyboard.Keyboard):
        comp.setAutoDraw(True)
    else:
        thisExp.addData('kb.started', comp.tStart)

def erase_comp(comp, t, tThisFlip, tThisFlipGlobal, frameN):
    """ Erases the stimuli displayed in trials. """
    comp.tStop = tThisFlip
    comp.tStopRefresh = tThisFlipGlobal
    comp.frameNStop = frameN
    comp.status = FINISHED
    if not isinstance(comp, keyboard.Keyboard):
        comp.setAutoDraw(False)
    else:
        thisExp.addData('kb.stopped', comp.tStop)

import random
import itertools

def generate_all_visuals(feature_config, target_pokemon, target_color):
    """
    Generate all visual stimuli for a feature condition run, including RSVP Pokémon sequences and peripheral stimulus grids.

    Args:
        feature_config (dict): Experiment parameters.
        target_pokemon (str): The Pokémon that will be the RSVP target.
        target_color (str): The color that will be the grid target.

    Returns:
        all_run_visuals (dict): For each run, a list of trial dicts with:
            'block_num', 'trial_num', 'vf', 'presentation_cond', 
            'periph_stim_grid', 'rsvp_sequence'.
        blanks_rsvps (dict): For each run, the blank block RSVP sequences.
    """
    # --- Extract config
    num_blocks = feature_config['NUM_BLOCKS']
    trials_per_block = feature_config['NUM_TRIALS']
    block_structure = feature_config['BLOCK_DESIGN']
    runs_total = feature_config['NUM_RUNS']
    blanks_per_run = feature_config['NUM_BLANKS']
    attention_conditions = feature_config['ATTENTION_CONDS']
    num_attention_types = len(attention_conditions)
    num_unique_runs = runs_total // num_attention_types # run visuals will be repeated once across each attention condition
    trials_per_run = num_blocks * trials_per_block
    run_duration = ((num_blocks/2)*trials_per_block)*SIM_TRIAL_DURATION + ((num_blocks/2)*trials_per_block)*SEQ_TRIAL_DURATION

    # --- Blank RSVP Generation ---
    def create_blank_rsvp_sequences():
        """Generate the RSVP lists for blank blocks."""
        num_unique_blanks = int(blanks_per_run * runs_total // 2)
        pokemon_per_blank = int(BLANK_BLOCK_DURATION // RSVP_RATE)
        all_blank_sequences = []
        for block in range(num_unique_blanks):
            # Get times for target occurrences in this block
            target_indices = []
            time_to_next_target = random.uniform(*POKEMON_TARGET_FREQ)
            while time_to_next_target < BLANK_BLOCK_DURATION:
                idx = int(time_to_next_target // RSVP_RATE)
                if idx < pokemon_per_blank:
                    target_indices.append(idx)
                time_to_next_target += random.uniform(*POKEMON_TARGET_FREQ)
            # Build sequence, prevent back-to-back repeats
            sequence = []
            for idx in range(pokemon_per_blank):
                if idx in target_indices:
                    sequence.append(target_pokemon)
                else:
                    distractor_options = [p for p in pokemon_names if p != target_pokemon]
                    if idx > 0:
                        distractor_options = [p for p in distractor_options if p != sequence[-1]]
                    sequence.append(random.choice(distractor_options))
            all_blank_sequences.append(sequence)
        # Repeat the sequences once per attention condition
        all_blank_sequences *= num_attention_types
        # Assign pairs (start/end) per run
        blank_rsvps = {}
        seq_idx = 0
        for run in range(1, runs_total + 1):
            blank_rsvps[run] = [all_blank_sequences[seq_idx], all_blank_sequences[seq_idx + 1]]
            seq_idx += num_attention_types
        return blank_rsvps

    # --- Peripheral Grid Generation ---
    def get_grid_assignments(target_color):
        """Return two lists: all grids with and without the target color."""
        all_color_combos = list(itertools.combinations(PERIPHERAL_STIM_COLORS, NUM_PSTIMS)) # all possible subsets of 4 colors from 6
        all_grids = []
        for combo in all_color_combos:
            all_grids.extend(itertools.permutations(combo)) # all possible ways to uniquely order those 4 colors
        random.shuffle(all_grids)
        grids_with_target = [grid for grid in all_grids if target_color in grid]
        grids_without_target = [grid for grid in all_grids if target_color not in grid]
        return grids_with_target, grids_without_target

    # --- Per-run trial dictionary generation ---
    def create_trials_for_run(run_id, grids_with_target, grids_without_target):
        """Create all trial dicts for a single run."""
        # Schedule trials that will contain the target color in their grid
        target_grid_trials = []
        first_target_trial = random.randint(*PSTIM_TARGET_FREQ)
        trial_idx = first_target_trial
        while trial_idx < trials_per_run:
            target_grid_trials.append(trial_idx)
            trial_idx += random.randint(*PSTIM_TARGET_FREQ)
        # Assign grids to each trial (unique if possible)
        grid_assignments = [None] * trials_per_run
        used_grids = set()
        for idx in range(trials_per_run):
            if idx in target_grid_trials:
                available = [g for g in grids_with_target if g not in used_grids] or grids_with_target
            else:
                available = [g for g in grids_without_target if g not in used_grids] or grids_without_target
            chosen_grid = random.choice(available)
            grid_assignments[idx] = chosen_grid
            used_grids.add(chosen_grid)
        # --- Pokémon RSVP Scheduling ---
        # When does the target Pokémon appear in the run?
        target_pokemon_times = []
        next_target_time = random.uniform(*POKEMON_TARGET_FREQ)
        while next_target_time < run_duration:
            target_pokemon_times.append(next_target_time)
            next_target_time += random.uniform(*POKEMON_TARGET_FREQ)
        # Map each target Pokémon onset to a trial & its time within trial
        trial_target_onsets = {i: [] for i in range(trials_per_run)}
        for t in target_pokemon_times:
            trial_num = int(t // SIM_TRIAL_DURATION)
            if trial_num < trials_per_run:
                time_within_trial = t - (trial_num * SIM_TRIAL_DURATION)
                trial_target_onsets[trial_num].append(time_within_trial)
        # --- Build trial dicts ---
        trial_dicts = []
        for block_idx, (visual_field, present_cond) in enumerate(block_structure):
            for trial_in_block in range(trials_per_block):
                trial_num = block_idx * trials_per_block + trial_in_block
                grid = grid_assignments[trial_num]
                # Build this trial's RSVP sequence
                num_pokemons_in_trial = int(SIM_TRIAL_DURATION // RSVP_RATE)
                sequence = []
                for step in range(num_pokemons_in_trial):
                    timepoint = step * RSVP_RATE
                    # Insert target Pokémon at scheduled times
                    if any(abs(timepoint - t) < RSVP_RATE / 2 for t in trial_target_onsets[trial_num]):
                        sequence.append(target_pokemon)
                    else:
                        options = [p for p in pokemon_names if p != target_pokemon]
                        if step > 0:
                            options = [p for p in options if p != sequence[-1]]
                        sequence.append(random.choice(options))
                trial_dicts.append({
                    'block_num': block_idx + 1,
                    'trial_num': trial_num,
                    'visual_field': visual_field,
                    'presentation_cond': present_cond,
                    'peripheral_grid': grid,
                    'rsvp_sequence': sequence
                })
        return trial_dicts

    # --- Main routine ---
    blanks_rsvps = create_blank_rsvp_sequences()
    grids_with_target, grids_without_target = get_grid_assignments(target_color)
    all_run_visuals = {}
    unique_run_trials = [create_trials_for_run(i + 1, grids_with_target, grids_without_target) for i in range(num_unique_runs)]
    for run_idx in range(1, runs_total + 1):
        # Use the same visuals for corresponding runs across attention conditions
        all_run_visuals[f"run{run_idx}"] = unique_run_trials[(run_idx - 1) % num_unique_runs]

    return all_run_visuals, blanks_rsvps
    
def run_trial(trial_dict, attention_cond, target_pokemon, target_color):
    """
    Function to run a single trial.

    Args:
       trial_dict(dict): Contains 'block_num', 'trial_num', 'visual_field', 'presentation_cond', 'peripheral_grid','rsvp_sequence'
       attention_cond (str): 'FIX' (looking for a pokemon) or 'COV' (looking for a colored circle)
       target_pokemon (str): Name of pokemon to look for
       target_color (str): Name of color to look for
    """

    # Reset all variables
    trialClock.reset()
    frameN = -1
    rsvp_index = 0
    pstim_index = 0
    target_pstim_onset = None
    target_pokemon_onset = None
    continueRoutine = True
    kb.clearEvents()
    kb.rt = []
    kb.keys = []
    kb_allKeys = []
    pstim_to_draw = []
    
    # Get trial parameters from the trial dictionary
    vf = trial_dict['visual_field']
    is_seq_trial = trial_dict['presentation_cond'] == 'SEQ'
    is_sim_trial = trial_dict['presentation_cond'] == 'SIM'
    pstim_grid = trial_dict['peripheral_grid']
    rsvp_sequence = trial_dict['rsvp_sequence']

    # Extract peripheral stim locations
    if vf[0] == 'R':
        grid_positions = [rvf_topright, rvf_topleft, rvf_botleft, rvf_botright]
    elif vf[0] == 'L':
        grid_positions = [lvf_topright, lvf_topleft, lvf_botleft, lvf_botright]

    # Map peripheral stim colors to their positions based on their order in pstim_grid
    for color, pos in zip(pstim_grid, grid_positions):
        pstim = visual.Circle(win, name = color, pos = pos, radius = PERIPHERAL_STIM_SIZE/2 , units = 'deg', anchor = 'center', fillColor=color, lineColor=color)
        pstim_dict = {color: pstim}
        pstim_to_draw.append(pstim) # add all the pstims to be drawn to a list
    
    # SEQ only: calculate the pstim onsets with ISI and add extra pokemon to the RSVP
    if is_seq_trial:
        pstim_onsets = [i*(PERIPH_STIM_DURATION+ISI) for i in range(len(pstim_to_draw))]
        additional_pokemon = math.ceil(NUM_PSTIMS * ISI / RSVP_RATE) 
        distractors = [pokemon_name for pokemon_name in pokemon_names if pokemon_name != target_pokemon] # only non-target pokemon will be added to the RSVP
        rsvp_sequence.extend(random.choices(distractors, k=additional_pokemon))
    #SIM only: select random start time for peripheral stimuli
    else:
        pstim_start_time = random.choice([0, 1, 2, 3])
    
    # Reset visual components before trial loop
    components = [pstim_dict, pokemon_dict, kb]
    for comp in components:
        if isinstance(comp, dict):
            for i in comp:
                comp[i].tStart = None
                comp[i].tStop = None
                comp[i].tStartRefresh = None
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
                
    # Save trial variables to the data file
    thisExp.addData('pres_cond', trial_dict['presentation_cond'])
    thisExp.addData('block', trial_dict['block_num'] )
    thisExp.addData('trial', trial_dict['trial_num']+1)
    thisExp.addData('visual_field', vf[0])
    thisExp.addData('peripheral_grid', pstim_grid)
    thisExp.addData('rsvp_stim', rsvp_sequence)

    # Start the trial loop
    while continueRoutine:
        t = trialClock.getTime() # t = 0 is start of the trial
        tThisFlip = win.getFutureFlipTime(clock=trialClock) # when the next screen flip will be, relative to the start of the trial
        tThisFlipGlobal = win.getFutureFlipTime(clock=globalClock) # when the next screen flip will be, relative to the start of the experiment
        frameN += 1

        if rsvp_index < len(rsvp_sequence): # make sure we don't go past the last pokemon in the rsvp
            current_pokemon = pokemon_dict[rsvp_sequence[rsvp_index]]
            is_final_pokemon = rsvp_index == len(rsvp_sequence) - 1 # check if this is the last pokemon in the RSVP

            # Visual flow for SIM condition
            if is_sim_trial:
                # Draw pokemon RSVP stream
                if current_pokemon.status == NOT_STARTED and tThisFlip >= 0: # Start RSVP stream at the start of the trial
                    draw_comp(current_pokemon, t, tThisFlip, tThisFlipGlobal, frameN)
                    current_pokemon_onset = current_pokemon.tStartRefresh # Record time when the pokemon was first presented (from global clock!)
                    if rsvp_sequence[rsvp_index] == target_pokemon and target_pokemon_onset == None:
                        thisExp.addData(f'{target_pokemon}.started', comp.tStart) # Add time when target pokemon appears (in reference to start of trial) to data file
                        target_pokemon_onset = current_pokemon_onset # If pokemon is target pokemon, save first onset time

                if current_pokemon.status == STARTED and tThisFlipGlobal > current_pokemon_onset + RSVP_RATE: # erase the current pokemon after rsvp rate elapses
                    erase_comp(current_pokemon, t, tThisFlip, tThisFlipGlobal, frameN)
                    if not is_final_pokemon:
                        rsvp_index += 1 # move to the next pokemon in the sequence unless this is the last one
                        pokemon_dict[rsvp_sequence[rsvp_index]].status = NOT_STARTED # reset status for the next pokemon to be drawn in case it has already been shown
                
                # Draw the peripheral stimulus grid all at once at a random timepoint
                for pstim in pstim_to_draw:
                    if pstim.status == NOT_STARTED and tThisFlip >= pstim_start_time:
                        draw_comp(pstim, t, tThisFlip, tThisFlipGlobal, frameN)
                        if pstim.color == target_color and target_pstim_onset == None:
                            target_pstim_onset = pstim.tStartRefresh # If target color in grid, save onset time (from global clock!)
                            thisExp.addData(f'{target_color}.started', comp.tStart) # Save target onset to data file (reference to start of trial)

                    if pstim.status == STARTED and tThisFlipGlobal > pstim.tStartRefresh + PERIPH_STIM_DURATION: # Erase grid after pstim duration
                        erase_comp(pstim, t, tThisFlip, tThisFlipGlobal, frameN)

            # Visual flow for SEQ condition
            if is_seq_trial:
                # Draw pokemon RSVP stream
                pokemon_duration = SEQ_TRIAL_DURATION - SIM_TRIAL_DURATION if is_final_pokemon else RSVP_RATE # last pokemon is displayed for 99ms in SEQ to account for ISI

                if current_pokemon.status == NOT_STARTED and tThisFlip >= 0: # Start RSVP stream at the start of the trial
                    draw_comp(current_pokemon, t, tThisFlip, tThisFlipGlobal, frameN)
                    current_pokemon_onset = current_pokemon.tStartRefresh # Record time when the pokemon was first presented
                    if rsvp_sequence[rsvp_index] == target_pokemon and target_pokemon_onset == None:
                        thisExp.addData(f'{target_pokemon}.started', comp.tStart) # Add time when target pokemon appears (in reference to start of trial) to data file
                        target_pokemon_onset = current_pokemon_onset # If pokemon is target pokemon, save first onset time

                if current_pokemon.status == STARTED and tThisFlipGlobal > current_pokemon_onset + pokemon_duration: # erase the current pokemon after rsvp rate elapses
                    erase_comp(current_pokemon, t, tThisFlip, tThisFlipGlobal, frameN)
                    if not is_final_pokemon:
                        rsvp_index += 1 # move to the next pokemon in the sequence unless this is the last one
                        pokemon_dict[rsvp_sequence[rsvp_index]].status = NOT_STARTED # reset status for the next pokemon to be drawn

                # Draw peripheral stimuli one at a time
                if pstim_index < len(pstim_to_draw):
                    current_pstim = pstim_to_draw[pstim_index]
                    scheduled_onset = pstim_onsets[pstim_index]
                    if current_pstim.status == NOT_STARTED and tThisFlip >= scheduled_onset:
                        draw_comp(current_pstim, t, tThisFlip, tThisFlipGlobal, frameN)
                        if current_pstim.color == target_color and target_pstim_onset == None:
                            target_pstim_onset = current_pstim.tStartRefresh # If target color shown, save onset time
                            thisExp.addData(f'{target_color}.started', comp.tStart) # Save target onset to data file (reference to start of trial)

                    if current_pstim.status == STARTED and tThisFlipGlobal > current_pstim.tStartRefresh + PERIPH_STIM_DURATION:
                        erase_comp(current_pstim, t, tThisFlip, tThisFlipGlobal, frameN)
                        pstim_index += 1

        # Start checking for key presses at start of trial
        if kb.status == NOT_STARTED and tThisFlip >= 0:
            draw_comp(kb, t, tThisFlip, tThisFlipGlobal, frameN)

        # Check keypresses once per frame
        elif kb.status == STARTED:
            keys = kb.getKeys(keyList=[RESPONSE_KEY, 'escape'], waitRelease=False)
            for key in keys:
                kb_allKeys.append(key)
                if key.name == 'escape':
                    end_task()
                if key.name == RESPONSE_KEY:
                    if not kb.keys:
                        thisExp.timestampOnFlip(win, 'kb.pressed') # based on the global clock
                        kb.keys = key.name
                        kb.rt = kb_allKeys[0].rt # based on start of trial (when kb is initialized)

        # End the trial after the total duration
        if (is_sim_trial and tThisFlip >= SIM_TRIAL_DURATION) or (is_seq_trial and tThisFlip >= SEQ_TRIAL_DURATION):
            continueRoutine = False 
            for comp in components: # make sure all components are erased
                if isinstance(comp, dict):
                    for i in comp:
                        erase_comp(comp[i], t, tThisFlip, tThisFlipGlobal, frameN)
                else:
                    erase_comp(comp, t, tThisFlip, tThisFlipGlobal, frameN)
            win.clearAutoDraw()
                    
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
        
        # Blank block at start of run; display an RSVP of distractors
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
        
        # Blank block at end of run; display an RSVP of distractors 
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

# Clear the window
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

# Start trial-level clock and clear the window
trialClock = core.Clock()
win.flip() 

# Run each feature condition
run_feature_cond(color_condition, 'Pikachu', 'red')

###### END EXPERIMENT #######################################################################################################

# Draw thank you text and run function to save data
thanks_text.draw()
win.flip()
keys = event.waitKeys(keyList=['space'])

end_task() 