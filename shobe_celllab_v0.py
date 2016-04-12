# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 08:39:05 2016

@author: Charlie Shobe

Modeling class cell lab code: grains dumped into a river

Conditions:
-Flow is from right to left
-All boundaries closed except right edge
-New grains are dumped in at the left edge each timestep
-grains are pulled by gravity, random walk due to turbulence, and are
moved left to right by flow (see transition probabilities)

This will produce a continuously updating plot of the river (in cross-section) and grains.
Once the model is finished, if you close figure 1, three plots of concentration
profiles at x = 50, x = 100, and x = 150 will be shown.
"""

import time
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from landlab import RasterModelGrid
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.oriented_raster_cts import OrientedRasterCTS

def setup_transition_list():
    """
    Creates and returns a list of Transition() objects to represent state
    transitions for a biased random walk, in which the rate of downward
    motion is greater than the rate in the other three directions.
    
    Parameters
    ----------
    (none)
    
    Returns
    -------
    xn_list : list of Transition objects
        List of objects that encode information about the link-state transitions.
    
    Notes
    -----
    State 0 represents fluid and state 1 represents a particle (such as a 
    sediment grain, tea leaf, or dissolved heavy particle).
    
    Try to add 3rd state, which is "bedrock" grain?
    
    The states and transitions are as follows:

    Pair state      Transition to       Process             Rate (cells/s)
    ==========      =============       =======             ==============
    0 (0-0)         (none)              -                   -
    1 (0-1)         2 (1-0)             left motion         0.1 #changed by CMS
    2 (1-0)         1 (0-1)             right motion        9.9
    3 (1-1)         (none)              -                   -
    4 (0-0)         (none)              -                   -
    5 (0-1)         2 (1-0)             down motion         10.55 #changed by cms
    6 (1-0)         1 (0-1)             up motion            9.45 #changed by cms
    7 (1-1)         (none)              -                   -
    """
    
    # Create an empty transition list
    xn_list = []
    
    # Append four transitions to the list.
    # Note that the arguments to the Transition() object constructor are:
    #  - Tuple representing starting pair state
    #    (left cell, right cell, orientation [0=horizontal])
    #  - Tuple representing new pair state
    #    (bottom cell, top cell, orientation [1=vertical])
    #  - Transition rate (cells per time step, in this case 1 sec)
    #  - Name for transition
    xn_list.append( Transition((0,1,0), (1,0,0), 0.1, 'left motion') )
    xn_list.append( Transition((1,0,0), (0,1,0), 9.9, 'right motion') )
    xn_list.append( Transition((0,1,1), (1,0,1), 10.55, 'down motion') )
    xn_list.append( Transition((1,0,1), (0,1,1), 9.45, 'up motion') )
    
    return xn_list
    
# INITIALIZE

# User-defined parameters
nr = 50  # number of rows in grid
nc = 200  # number of columns in grid
plot_interval = 1.0   # time interval for plotting, sec #THIS IS MODEL TIMESTEP
run_duration = 20.0   # duration of run, sec
report_interval = 10.0  # report interval, in real-time seconds

# Remember the clock time, and calculate when we next want to report
# progress.
current_real_time = time.time()
next_report = current_real_time + report_interval

# Create grid
mg = RasterModelGrid(nr, nc, 1.0)

# Make the boundaries be walls #right, bottom, left, top
mg.set_closed_boundaries_at_grid_edges(False, True, True, True)

# Set up the states and pair transitions.
ns_dict = { 0 : 'fluid', 1 : 'particle' }
xn_list = setup_transition_list()

# Create the node-state array and attach it to the grid
node_state_grid = mg.add_zeros('node', 'node_state_map', dtype=int)

# Initialize the node-state array: here, the initial condition is a pile of
# resting grains at the bottom of a container.
left_cols = np.where(mg.node_x<0.05*nc)[0]
node_state_grid[left_cols] = 1

# For visual display purposes, set all boundary nodes to fluid
node_state_grid[mg.closed_boundary_nodes] = 0

# Create the CA model
ca = OrientedRasterCTS(mg, ns_dict, xn_list, node_state_grid)

grain = '#5F594D'
fluid = '#D0E4F2'
clist = [fluid,grain]
my_cmap = matplotlib.colors.ListedColormap(clist)

# Create a CAPlotter object for handling screen display
ca_plotter = CAPlotter(ca, cmap=my_cmap)

# Plot the initial grid
ca_plotter.update_plot()

# RUN
current_time = 0.0
while current_time < run_duration:
    
    # Once in a while, print out simulation and real time to let the user
    # know that the sim is running ok
    current_real_time = time.time()
    #CMS says: supply new grains at the left side
    node_state_grid[left_cols] = 1   
    if current_real_time >= next_report:
        print 'Current sim time',current_time,'(',100*current_time/run_duration,'%)'
        next_report = current_real_time + report_interval
    
    # Run the model forward in time until the next output step
    ca.run(current_time+plot_interval, ca.node_state, 
           plot_each_transition=False)
    current_time += plot_interval
    
    # Plot the current grid
    ca_plotter.update_plot()

# FINALIZE

# Plot
print 'Model complete; close figure for concentration profiles'
ca_plotter.finalize()
plt.close(1)
#plt.savefig('final_grain_config.png')

# Calculate concentration profiles
c50 = np.zeros(nr)
left_bounds = np.arange(50 - 10, nr * nc, nc)
right_bounds = np.arange(50 + 10, nr * nc, nc)
for r in range(nr):
    c50[r] = np.mean(node_state_grid[left_bounds[r]:right_bounds[r]])
plt.figure(2)
plt.plot(c50, range(nr), 'o')
plt.xlabel('Concentration [cells / 20 cell window]')
plt.ylabel('Height [cells]')  
plt.xlim(0, 1)
plt.title('Concentration Profile at Position 50')

c100 = np.zeros(nr)
left_bounds = np.arange(100 - 10, nr * nc, nc)
right_bounds = np.arange(100 + 10, nr * nc, nc)
for r in range(nr):
    c100[r] = np.mean(node_state_grid[left_bounds[r]:right_bounds[r]])
plt.figure(3)
plt.plot(c100, range(nr), 'o')
plt.xlabel('Concentration [cells / 20 cell window]')
plt.ylabel('Height [cells]')    
plt.xlim(0, 1)
plt.title('Concentration Profile at Position 100')

c150 = np.zeros(nr)
left_bounds = np.arange(150 - 10, nr * nc, nc)
right_bounds = np.arange(150 + 10, nr * nc, nc)
for r in range(nr):
    c150[r] = np.mean(node_state_grid[left_bounds[r]:right_bounds[r]])
plt.figure(4)
plt.plot(c150, range(nr), 'o')
plt.xlabel('Concentration [cells / 20 cell window]')
plt.ylabel('Height [cells]')  
plt.xlim(0, 1)
plt.title('Concentration Profile at Position 150')
plt.show()
