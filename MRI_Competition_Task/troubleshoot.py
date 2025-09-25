import random
import itertools

# Experiment design
ATTENTION_CONDS = ['FIX', 'COV']
NUM_RUNS = 6 # per feature condition
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

pokemon_names = ["Bulbasaur", "Pikachu", "Squirtle", "Charmander", "Magikarp", "Raticate", "Pidgey",
    "Metapod", "Jigglypuff", "Butterfree", "Psyduck", "Caterpie", "Krabby",
    "Haunter", "Vulpix", "Eevee", "Sandshrew", "Wartortle", "Charmeleon", "Clefairy",
    "Ponyta", "Mankey"]

target_pokemon = 'Pikachu'
target_color = 'red'

def generate_blank_rsvps():
    """" Generate RSVP sequences for all blank blocks.
    Returns a dict where the keys are the unique run indices (0-indexed) and the values are RSVP sequences for pre and post blank blocks."""
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
    Returns a dictionary where the keys are the unique run indices (0-indexed) and the values are lists 
    of each trial's RSVP sequence for that run.
        
    Example:
        Calling generate_trial_rsvps()[0][0] will give the RSVP sequence (list of pokemon names) 
        of the first trial of the first run. 
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
    """ For each run, assign each trial a peripheral stimulus grid. Returns a list of lists (one list of grids per run)."""

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
    Create a dictionary of the parameters needed to run a trial.

    Args:
        run_idx (int): index of the run (0-indexed)
        all_grids (list): all of the grid assignments for all unique runs.
        all_trial_rsvps (dict): all of the rsvp sequences for all trials for all unique runs. 

    Returns:
        trial_dicts (list): Each value is a trial dictionary including:
            'block_num', 'trial_num', 'visual_field', 'presentation_cond', 
            'peripheral_grid', 'rsvp_seq'. 
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

all_grids = assign_grids()
all_trial_rsvps = generate_trial_rsvps()

run1 = create_trial_dicts(0, all_grids, all_trial_rsvps)

import pandas as pd

pokemon = [
    'Squirtle', 'Sandshrew', 'Charmander', 'Metapod', 'Clefairy', 'Eevee',
    'Squirtle', 'Haunter', 'Magikarp', 'Squirtle', 'Charmeleon', 'Sandshrew',
    'Mankey', 'Bulbasaur', 'Butterfree', 'Pidgey', 'Magikarp', 'Ponyta',
    'Butterfree', 'Charmander', 'Eevee', 'Charmeleon', 'Charmander', 'Vulpix',
    'Pikachu', 'Charmander', 'Raticate', 'Charmander', 'Sandshrew', 'Krabby',
    'Squirtle', 'Bulbasaur', 'Clefairy', 'Krabby', 'Metapod', 'Krabby',
    'Ponyta', 'Sandshrew', 'Haunter', 'Ponyta', 'Krabby', 'Squirtle',
    'Pidgey', 'Mankey', 'Jigglypuff', 'Ponyta', 'Sandshrew', 'Charmeleon',
    'Vulpix', 'Jigglypuff', 'Pikachu', 'Clefairy', 'Charmander', 'Eevee',
    'Bulbasaur', 'Vulpix', 'Wartortle', 'Vulpix', 'Pidgey', 'Clefairy',
    'Psyduck', 'Eevee', 'Mankey', 'Sandshrew'
]

# Create a single-column DataFrame
df = pd.DataFrame(pokemon, columns=['Pokemon'])

# Save to Excel
df.to_excel("pokemon_list.xlsx", index=False)
