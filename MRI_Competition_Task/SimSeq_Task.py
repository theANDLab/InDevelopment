from psychopy import prefs, plugins, sound, gui, visual, core, data, event, logging, clock, colors, layout, hardware, monitors
from psychopy.constants import (NOT_STARTED, STARTED, STOPPED, FINISHED, PRESSED, RELEASED, FOREVER, priority)
plugins.activatePlugins()
prefs.hardware['audioLib'] = 'ptb'
prefs.hardware['audioLatencyMode'] = '3'
from psychopy.hardware import keyboard
import math
import numpy as np
import csv
import itertools
import random
import itertools
import os
import sys
import time
logging.console.setLevel(logging.ERROR)

###### EXPERIMENT PARAMETERS #######################################################################################################

# Initialize the global clock and keyboard
globalClock = core.Clock()
kb = keyboard.Keyboard(clock = globalClock)

# Experiment design
ATTENTION_CONDS = ['FIX', 'COV']
NUM_RUNS = 6 # number of runs per feature condition
NUM_SIM_BLOCKS = 6 # per run
NUM_SEQ_BLOCKS = 6 # per run
NUM_BLANK_BLOCKS = 2 # per run
NUM_TRIALS = 3 # per block
BLOCK_DESIGN = [('RVF','SIM'),('LVF','SEQ'),('RVF','SEQ'),('LVF','SIM'),('RVF','SEQ'),('LVF','SIM'),
                ('RVF','SIM'),('LVF','SEQ'),('RVF','SIM'),('LVF','SEQ'),('RVF','SEQ'),('LVF','SIM')]

# Size
PERIPHERAL_STIM_SIZE = 1.75 #DVA; size of each peripheral stimulus (circle)
POKEMON_SIZE = [1.5, 1.5] # DVA, size of the pokemon during RSVP
POKEMON_POS = (0,0) # location of rsvp pokemon
NUM_PSTIMS = 4 # number of peripheral stimuli
GRID_SIZE = 4 #DVA; height and width of the peripheral stimulus grid
ECCENTRICITY = 7 #DVA from the center of the grid to the center of each peripheral stimulus
PERIPHERAL_STIM_COLORS = ['red', 'blue', 'green', 'yellow', 'magenta', 'cyan'] 

# Timing
BLANK_BLOCK_DURATION = 16 # seconds
PERIPH_STIM_DURATION = 1 # seconds
RSVP_RATE = 0.25 #sec; duration of RSVP pokemon presentation
TRIAL_DURATION = PERIPH_STIM_DURATION*NUM_PSTIMS # sec
PSTIM_TARGET_FREQ = [1,3] # pstim color targets will occur every 1-3 trials
POKEMON_TARGET_FREQ = [15,30] # pokemon targets will occur every 15-30 pokemon (3.75-7.5s) in the RSVP
RESPONSE_WINDOW = 1.5 # sec; responses after this will be coded as FAs

# Response key
RESPONSE_KEY = 'space'

###### EXPERIMENT SETUP #######################################################################################################

exp_name = 'SimSeq-DV'
exp_info = {'Participant ID': '9999', 
            'Session': '001',
            'Pokemon': 'Pikachu',
            'Color': 'red'
            }
dlg = gui.DlgFromDict(dictionary=exp_info, title=exp_name)
if dlg.OK == False:
    core.quit()  

# Get targets
target_pokemon = exp_info['Pokemon'].strip().capitalize() # ensures first letter is capitalized
target_color = exp_info['Color']

# Establish data output directory
time_str = time.strftime("_%m_%d_%Y_%H_%M", time.localtime())
root_dir = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(root_dir, 'data', f"{exp_name}_{exp_info['Participant ID']}_Session{exp_info['Session']}_{time_str}")
os.makedirs(data_folder, exist_ok=True)
filename = os.path.join(data_folder, f"{exp_name}_{exp_info['Participant ID']}_Session{exp_info['Session']}")
blankblock_filename = os.path.join(data_folder, f"blank_blocks_{exp_info['Participant ID']}_Session{exp_info['Session']}")

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

blankExp = data.ExperimentHandler(name='blank_blocks',extraInfo=exp_info,
                                savePickle=True, saveWideText=True,
                                dataFileName=blankblock_filename)

###### INITIALIZE VISUAL COMPONENTS #######################################################################################################

# Load Pokémon images from folder
pokemon_dir = os.path.join(root_dir, 'pokemon_lightgray')
pokemon_names = ["Bulbasaur", "Pikachu", "Squirtle", "Charmander", "Magikarp", "Raticate", "Pidgey",
    "Metapod", "Jigglypuff", "Butterfree", "Psyduck", "Caterpie", "Krabby",
    "Haunter", "Vulpix", "Eevee", "Sandshrew", "Wartortle", "Charmeleon", "Clefairy",
    "Ponyta", "Mankey"]
    
pokemon_dict = {name: visual.ImageStim(win, name=name, image=os.path.join(pokemon_dir, f"{i+1:03}.png"))
    for i, name in enumerate(pokemon_names)} # dictionary of the pokemon where the key is their name and the values are the ImageStim

# Text components for welcome, instructions, and end screens
welcome_text = visual.TextStim(win, pos=(0,0), height= 1.5, wrapWidth=27, text=("Welcome to the Pokémon Party game!"))
fix_instructions_text = visual.TextStim(win, pos = (0,0), wrapWidth=27, text=(
        "There's a Pokémon Party happening right now, and the Pokémon are playing hide and seek!\n\n"
        f"The Pokémon are having trouble finding {target_pokemon}! Can you help them?\n\n"
        f"Press the button as fast as you can every time you see {target_pokemon}.\n\n\n"
        "Ready to start playing?"
    ))
cov_instructions_text = visual.TextStim(win, pos = (0,0), wrapWidth=27, text=(
        "There's a Pokémon Party happening right now, and the Pokémon are getting hungry!\n\n"
        f"The Pokémon like to eat {target_color} circles! Can you help feed them?\n\n"
        f"Press the button as fast as you can every time you see a {target_color} circle.\n\n\n"
        "Ready to start playing?"
    ))
end_text = visual.TextStim(win, wrapWidth=27, text=())
thanks_text = visual.TextStim(win, wrapWidth=30, text=("Thanks for coming to the Pokémon Party!"))

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

def erase_comp(comp, t, tThisFlip, tThisFlipGlobal, frameN):
    """ Erases the stimuli displayed in trials. """
    comp.tStop = tThisFlip
    comp.tStopRefresh = tThisFlipGlobal
    comp.frameNStop = frameN
    comp.status = FINISHED
    if not isinstance(comp, keyboard.Keyboard):
        comp.setAutoDraw(False)
        
# TODO: function to save visuals
    
def generate_blank_rsvps():
    """" 
    Returns a dictionary where the keys are the unique run indices and the values are lists 
    of each blank block's RSVP sequence for that run. (0-indexed)
    
    Example:
        Calling generate_blank_rsvps()[2][0] will give the RSVP sequence (list of pokemon names) 
        of the first blank block in the third run. 
    """
    num_unique_runs = int(NUM_RUNS//len(ATTENTION_CONDS))
    num_unique_blanks = NUM_BLANK_BLOCKS * num_unique_runs
    pokemon_per_blank = int(BLANK_BLOCK_DURATION // RSVP_RATE)
    all_blank_sequences = []
    
    for blank_block in range(num_unique_blanks):
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
        all_blank_sequences.append(sequence)

    # Assign pairs (start/end) per run
    all_blank_rsvps = {}
    seq_idx = 0
    for run in range(num_unique_runs):
        all_blank_rsvps[run] = [all_blank_sequences[seq_idx], all_blank_sequences[seq_idx + 1]]
        seq_idx += len(ATTENTION_CONDS)
        
    return all_blank_rsvps
    
def generate_trial_rsvps():
    """ 
    Returns a dictionary where the keys are the unique run indices and the values are lists 
    of each trial's RSVP sequence for that run.(0-indexed)
        
    Example:
        Calling generate_trial_rsvps()[0][32] will give the RSVP sequence (list of pokemon names) 
        of the 33rd trial in the first run. 
    """

    pokemon_per_trial = int(TRIAL_DURATION // RSVP_RATE)
    total_pokemon = int((NUM_TRIALS*(NUM_SIM_BLOCKS+NUM_SEQ_BLOCKS))*pokemon_per_trial) # in one run
    run_seq_list = []

    for unique_run in range(int(NUM_RUNS//len(ATTENTION_CONDS))):
        target_pokemon_idx = []
        next_target_idx = random.randint(*POKEMON_TARGET_FREQ)
        while next_target_idx < total_pokemon:
            target_pokemon_idx.append(next_target_idx)
            next_target_idx += random.randint(*POKEMON_TARGET_FREQ) # list of all target for the run
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
        run_seq_list.append(run_sequence)

    all_trial_rsvps = {}
    for run_idx, run_sequence in enumerate(run_seq_list):
        trial_sequences = [run_sequence[i:i + pokemon_per_trial] for i in range(0, len(run_sequence), pokemon_per_trial)]
        all_trial_rsvps[run_idx] = trial_sequences

    return all_trial_rsvps
    
def assign_grids():
    """ 
    Returns a list of lists of each trial's RSVP sequence for that run.
        
    Example:
        Calling assign_grids()[0][32] will give the grid layout (list of colors) for 
        the 33rd trial in the first run. 
    """

    trials_per_run = NUM_TRIALS*(NUM_SEQ_BLOCKS+NUM_SIM_BLOCKS)
    num_attention_conds = len(ATTENTION_CONDS)

    # Get lists of target and nontarget grids
    all_color_combos = list(itertools.combinations(PERIPHERAL_STIM_COLORS, NUM_PSTIMS)) # all possible subsets of 4 colors from 6
    all_configs = []
    for combo in all_color_combos:
        all_configs.extend(itertools.permutations(combo)) # all possible ways to uniquely order those 4 colors
    random.shuffle(all_configs)
    rvf_target_grids = [grid for grid in all_configs if target_color in grid and grid[2] == target_color] # target in bottom left
    lvf_target_grids = [grid for grid in all_configs if target_color in grid and grid[3] == target_color] # target in bottom right
    grids_without_target = [grid for grid in all_configs if target_color not in grid]
    
    all_grids = []
    for run in range(NUM_RUNS//num_attention_conds):
        # Create list of trial indices for trials that will contain the target color
        target_trials = []
        target_trial_idx = random.randint(*PSTIM_TARGET_FREQ) # randomly select first target trial
        while target_trial_idx <= trials_per_run:
            target_trials.append(target_trial_idx)
            target_trial_idx += random.randint(*PSTIM_TARGET_FREQ)
            
        # Assign grids to each trial (unique if possible)
        grid_assignments = []
        used_grids = set()
        for trial_idx in range(trials_per_run):
            if trial_idx in target_trials:
                lvf_available = [g for g in lvf_target_grids if g not in used_grids] or lvf_target_grids
                rvf_available = [g for g in rvf_target_grids if g not in used_grids] or rvf_target_grids
            else:
                available = [g for g in grids_without_target if g not in used_grids] or grids_without_target
            # If target trial, assign grid based on the visual field of the block that the trial falls within
            if trial_idx in target_trials:
                block_num = (trial_idx// NUM_TRIALS) % len(BLOCK_DESIGN)
                if BLOCK_DESIGN[block_num][0] == 'RVF':
                    chosen_grid = random.choice(rvf_available)
                else:
                    chosen_grid = random.choice(lvf_available)
            else:
                chosen_grid = random.choice(available)
            grid_assignments.append(chosen_grid)
            used_grids.add(chosen_grid)
        all_grids.append(grid_assignments)

    return all_grids
    
def create_trial_dicts(run_idx, all_grids, all_trial_rsvps):
    """
    Generates trial dictionaries for all of the trials in a run. Call once at the start of each run. 

    Args:
        run_idx (int): index of the run (0-indexed)
        all_grids (list): all of the grid assignments for all unique runs.
        all_trial_rsvps (dict): all of the rsvp sequences for all trials for all unique runs. 
        
    Returns:
        trial_dicts (list): Each value is a trial dictionary including:
            'block_num', 'trial_num', 'visual_field', 'presentation_cond', 
            'peripheral_grid', 'rsvp_seq'. 
            
    Example:
        Calling create_trial_dicts(0, all_grids_all_trial_rsvps)[32] gives the trial dictionary of
        the 33rd trial in the first run. 
    """
    trial_dicts = []
    for block_idx, (visual_field, present_cond) in enumerate(BLOCK_DESIGN):
        for trial_in_block in range(NUM_TRIALS):
            trial_idx = block_idx * NUM_TRIALS + trial_in_block
            trial_dicts.append({
                'block_num': block_idx + 1,
                'trial_num': trial_idx + 1,
                'visual_field': visual_field,
                'presentation_cond': present_cond,
                'peripheral_grid': all_grids[run_idx][trial_idx],
                'rsvp_seq': all_trial_rsvps[run_idx][trial_idx]
            }) 
    return trial_dicts
    
def run_blank_block(rsvp, run_num, blank_num): 
    """ Function to run a single blank block. """
    
    # Add data to data file
    blankExp.addData('blank_block.started', globalClock.getTime(format='float'))
    blankExp.addData('rsvp_stim', rsvp)
    
    # Reset tracking variables
    last_target_onset = None # will store the most recent target onset
    total_hits = 0
    
    for current_pokemon in rsvp:
        blankExp.addData('run', run_num)
        blankExp.addData('blank_block', blank_num)
        hit = 0
        is_target = (current_pokemon == target_pokemon)
        
        # Set pokemon location
        pokemon_dict[current_pokemon].pos = (0,0)

        # Get a time stamp of when the while loop started
        start_time = globalClock.getTime(format='float')
        blankExp.addData('loop_start', start_time)
        pokemon_onset = None
        
        # while loop will draw pokemon and check for keypresses for RSVP_RATE duration
        while (globalClock.getTime() - start_time) < RSVP_RATE:
            
            # Draw pokemon every frame
            pokemon_dict[current_pokemon].draw()
            
            # Capture pokemon onset once during first window flip
            if pokemon_onset is None:
                win.callOnFlip(kb.clearEvents)
                def on_flip():
                    nonlocal pokemon_onset, last_target_onset
                    t = globalClock.getTime(format='float')
                    pokemon_onset = t
                    blankExp.addData('pokemon', current_pokemon)
                    blankExp.addData('onset', t)
                    if is_target:
                        last_target_onset = t
                        blankExp.addData('target_onset', t)
                win.callOnFlip(on_flip)
            win.flip()
            
            # Get all keys pressed during that 250ms window
            keys = kb.getKeys(keyList=[RESPONSE_KEY, 'escape'], waitRelease=False, clear = False)
            
            # Analyze key presses
            for key in keys:
                if key.name == 'escape':
                    end_task()
                    
                if key.name == RESPONSE_KEY:
                    press_time = key.rt # relative to globalClock, so NOT rt from target onset
                    blankExp.addData('pokemon', pokemon)
                    blankExp.addData('press_time', press_time)
                        
                    # Score against the most recent target onset
                    if last_target_onset is not None:
                        # compute rt since last target onset
                        rt = press_time - last_target_onset
                        blankExp.addData('rt', rt)
                        if rt <= RESPONSE_WINDOW:
                            hit = 1
                            total_hits += hit
                            blankExp.addData('hit', hit)
                    # next row in data file after a keypress to ensure we capture all of them
        blankExp.nextEntry()
    
    # Save block data
    blankExp.addData('blank_block.stopped', globalClock.getTime(format='float'))
    blankExp.nextEntry()
    
    return total_hits
    
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
    rsvp_sequence = trial_dict['rsvp_seq']

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

    pstim_onsets = [0, 1, 2, 3] # seconds; each peripheral stimuli will be shown one at a time in SEQ condition
    if is_sim_trial:
        pstim_start_time = random.choice(pstim_onsets) # select random start time for grid in sim trials

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
    thisExp.addData('trial', trial_dict['trial_num'])
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
            
            # Draw pokemon RSVP stream
            if current_pokemon.status == NOT_STARTED and tThisFlip >= 0: # Start RSVP stream at the start of the trial
                draw_comp(current_pokemon, t, tThisFlip, tThisFlipGlobal, frameN)
                current_pokemon_onset = current_pokemon.tStartRefresh # Record time when the pokemon was first presented
                if current_pokemon.name == target_pokemon and target_pokemon_onset == None:
                    thisExp.addData(f'{current_pokemon.name}.started', current_pokemon.tStart) # Add time when target pokemon appears (in reference to start of trial) to data file
                    thisExp.addData('target_onset2', current_pokemon_onset)
                    target_pokemon_onset = current_pokemon_onset # If pokemon is target pokemon, save first onset time

            if current_pokemon.status == STARTED and tThisFlipGlobal > current_pokemon_onset + RSVP_RATE: # erase the current pokemon after rsvp rate elapses
                erase_comp(current_pokemon, t, tThisFlip, tThisFlipGlobal, frameN)
                if not is_final_pokemon:
                    rsvp_index += 1 # move to the next pokemon in the sequence unless this is the last one
                    pokemon_dict[rsvp_sequence[rsvp_index]].status = NOT_STARTED # reset status for the next pokemon to be drawn
                   
            # SIM trials: draw the peripheral stimulus grid all at once at a random timepoint
            if is_sim_trial:                
                for pstim in pstim_to_draw:
                    if pstim.status == NOT_STARTED and tThisFlip >= pstim_start_time:
                        draw_comp(pstim, t, tThisFlip, tThisFlipGlobal, frameN)
                        if pstim.color == target_color and target_pstim_onset == None:
                            target_pstim_onset = pstim.tStartRefresh # If target color in grid, save onset time (from global clock!)
                            thisExp.addData(f'{target_color}.started', comp.tStart) # Save target onset to data file (reference to start of trial)
                    if pstim.status == STARTED and tThisFlipGlobal > pstim.tStartRefresh + PERIPH_STIM_DURATION: # Erase grid after pstim duration
                        erase_comp(pstim, t, tThisFlip, tThisFlipGlobal, frameN)
            
            # SEQ trials: Draw peripheral stimuli one at a time
            if is_seq_trial:
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

        # End the trial after the total duration
        if tThisFlip >= TRIAL_DURATION:
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
        elif first_press and first_press.rt >= target_onset: # key pressed after target was displayed
            rt = first_press.rt - target_onset
            print("target shown", "RT:", first_press.rt, "target onset time:", target_onset)
            if rt < RESPONSE_WINDOW: # pressed within response window -> hit
                accuracy = 1
                trial_outcome = "hit"
            else: # pressed after response window -> FA
                accuracy = 0
                trial_outcome = "false alarm"
        elif first_press and first_press.rt < target_onset: # key pressed before target displayed -> FA
            print("target shown", "RT:", first_press.rt, "target onset time:", target_onset)
            accuracy = 0
            trial_outcome = "false alarm"
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

def perform_one_run(feat_cond, run_idx, attention_cond, blanks_rsvps, all_grids, trial_rsvps, target_pokemon, target_color):
    
    # Save attention condition and run start time to data file
    thisExp.addData('attention_cond', attention_cond)
    thisExp.addData('run.started', globalClock.getTime(format='float'))
    
    # Reset run accuracy
    run_accuracy = 0
    
    # Show instructions for the attention condition
    if attention_cond == "FIX":
        fix_instructions_text.draw()
    if attention_cond == "COV":
        cov_instructions_text.draw()
    win.flip()
    
    # Move on from instructions when escape or space is pressed
    keys = event.waitKeys(keyList=['space', 'escape'])
    if 'escape' in keys:
        end_task()
    thisExp.addData('instructions.stopped', globalClock.getTime(format='float')) 
    
    # Extract visuals for this run
    trial_dicts = create_trial_dicts(run_idx, all_grids, trial_rsvps) # create trial dicts for this run
    blank_block_rsvps = blanks_rsvps[run_idx] # 2 RSVP lists for the blank blocks in this run
    
    # Run blank block before trials
    thisExp.addData('blank_block.started', globalClock.getTime(format='float'))
    run_blank_block(blank_block_rsvps[0], run_idx+1, 1)
    thisExp.addData('blank_block.stopped', globalClock.getTime(format='float'))

#    # Run all trials using trial dictionaries, updating accuracy
#    for trial_dict in trial_dicts:
#        run_accuracy += run_trial(trial_dict, attention_cond, target_pokemon, target_color)
    
    # Run blank block after trials
    thisExp.addData('blank_block.started', globalClock.getTime(format='float'))
    run_blank_block(blank_block_rsvps[1], run_idx+1, 2)
    thisExp.addData('blank_block.stopped', globalClock.getTime(format='float'))
    
    # Show feedback statement at the end of the run based on attention condition
    if attention_cond == 'FIX':
        end_text.text = (f"Great job finding {target_pokemon}!")
        pokemon_dict[target_pokemon].pos = (0, -5) 
        pokemon_dict[target_pokemon].size = (5,5)
        pokemon_dict[target_pokemon].draw()
    elif attention_cond == 'COV':
        end_text.text = (f"Great job feeding the Pokémon {target_color} circles!")
        target_pstim = visual.Circle(win, name = target_color, pos = (0,-5), radius = 5/2 , units = 'deg', anchor = 'center', fillColor=target_color, lineColor=target_color)
        target_pstim.draw()
    end_text.draw()
    win.flip()
    
    # Reset size and position of the target pokemon
    pokemon_dict[target_pokemon].size = POKEMON_SIZE
    pokemon_dict[target_pokemon].pos = (0,0)
    
    keys = event.waitKeys(keyList=['space', 'escape'])
    if 'escape' in keys:
        end_task()
 
def run_feature_cond(condition):
    """
    Runs all of the blocks in a single feature condition. 
    
    Args:
        feature_condition (str): color, motion, or color_motion
    """
    
    # Step 1: Generate the RSVP sequences and peripheral stimulus grids for all runs in the feature condition
    blanks_rsvps = generate_blank_rsvps()
    all_grids = assign_grids()
    trial_rsvps = generate_trial_rsvps()

    # Step 2: Perform 3 Color-FIX and 3 Color-COV runs
    for attention_cond in ATTENTION_CONDS:
        for run in range(NUM_RUNS//len(ATTENTION_CONDS)):
            perform_one_run(condition, run, attention_cond, blanks_rsvps, all_grids, trial_rsvps, target_pokemon, target_color)
    
###### WELCOME SCREEN #######################################################################################################

# Clear the window
win.flip() 

# Mark the experiment as started
exp_info['expStart'] = data.getDateStr(format = '%Y-%m-%d %Hh%M.%S.%f %z', fractionalSecondDigits=6)
thisExp.status = STARTED

# Draw welcome screen with Pokémon images
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

# Run color feature condition
run_feature_cond('color')

###### END EXPERIMENT #######################################################################################################

# Draw thank you text and run function to save data
thanks_text.draw()
win.flip()
keys = event.waitKeys(keyList=['space'])

end_task() 