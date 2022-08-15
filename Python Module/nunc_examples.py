# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 19:02:54 2022

@author: austine
"""
import numpy as np
import seqChange.np.NUNC as NUNC

'''
Contents
----------------

Example 1:
    Example of NUNC being used to detect a change, and then stopping further 
    testing.
Example 2:
    Demonstrates how increasing the test threshold prevents change detection.
Example 3:
    Shows an offline example where all points are tested and an indicator is
    returned for no change.
Example 4:
    Illlustrates NUNC on a data stream with multiple changes, continuing to 
    test after each change is declared.
Example 5: 
    Addresses the issue in Example 4 where the same changepoints are identified
    multiple times.
    

Introduction
------------------
NUNC can be used to monitor a sequence of data for a changepoint.
To do this we first instantiate the class nunc_local_state with a data point
and window size and then sequentially update the state with new observations.
We declare a change if the cost crosses a threshold.

The costs in the nunc_local_state are the costs in the window of data, so
to find the change we check if the max of the costs exceeds and threshold,
and then also use the location of that maximum. To do this we make use of the
function argmax (see below)

'''
    
def argmax(l) :
    #this function is needed to compute the location of the max when 
    #performing each sequential update
    pos = max(range(len(l)),key=lambda i: l[i])
    return (l[pos],pos) 

#In this first example we show that a data stream with a change being 
#monitored. A change is detected at point 52 at time 67.
#In this example we stop testing once a change is declared. For a different
#approach consider example 4.
#Note that there will be some uncertainty about the exact location
#of the change.

np.random.seed(1)
X = np.concatenate([np.random.normal(0, 1, 50), np.random.normal(5,5,150)])
X = np.ndarray.tolist(X)

w = 50 #window size
K = 3 #quantiles to use
S = NUNC.nunc_local_state(X[0], w) #initial state

threshold = 10
dtime = 1 #number of points observed
for x in X[1:]:
    S.update(K, x) #update state with new point
    dtime += 1
    if len(S.window) == w: #once the window is full start testing for changes
        (max_cost, position) = argmax(S.costs) #find max cost and position
        if max_cost > threshold: #if max cost exceeds threshold declare change
            print(dtime, dtime - w + position)
            break
    #above is the time of change and location in data stream
    
#In the second example we show how nothing is returned if the threshold
#is made larger.

np.random.seed(2)
X = np.concatenate([np.random.normal(0, 1, 50), np.random.normal(5,5,150)])
X = np.ndarray.tolist(X)

w = 50 #window size
K = 3 #quantiles to use
S = NUNC.nunc_local_state(X[0], w) #initial state

threshold = 30
dtime = 1 #number of points observed
for x in X[1:]:
    S.update(K, x) #update state with new point
    dtime += 1
    if len(S.window) == w: #once the window is full start testing for changes
        (max_cost, position) = argmax(S.costs) #find max cost and position
        if max_cost > threshold: #if max cost exceeds threshold declare change
            print(dtime, dtime - w + position)
            break
    #above is the time of change and location in data stream
    
#In this third example we consider how - in an offline setting - it may assist
#a user (when the data stream length is known) to record a lack of change 
#detection with a -1:
    
np.random.seed(3)
X = np.random.normal(0, 1, 200)
X = np.ndarray.tolist(X)

w = 50 #window size
K = 3 #quantiles to use
S = NUNC.nunc_local_state(X[0], w) #initial state

threshold = 30
dtime = 1 #number of points observed
for x in X[1:]:
    S.update(K, x) #update state with new point
    dtime += 1
    if len(S.window) == w: #once the window is full start testing for changes
        (max_cost, position) = argmax(S.costs) #find max cost and position
        if max_cost > threshold: #if max cost exceeds threshold declare change
            print(dtime, dtime - w + position)
            break
    print(dtime, -1) #no change detected
    #above is the time of change and location in data stream
    
#In this fourth example we show how NUNC can keep detecting changes rather 
#than breaking after the first one
    
np.random.seed(1)
X = np.concatenate([np.random.normal(0, 1, 50), np.random.normal(5,5,150),
                   np.random.normal(20, 1, 100), np.random.poisson(3, 100)])
X = np.ndarray.tolist(X)

w = 50 #window size
K = 3 #quantiles to use
S = NUNC.nunc_local_state(X[0], w) #initial state
change_locations = [] #for storing changepoints and detection times
detection_times = []

threshold = 10
dtime = 1 #number of points observed
for x in X[1:]:
    S.update(K, x) #update state with new point
    dtime += 1
    if len(S.window) == w: #once the window is full start testing for changes
        (max_cost, position) = argmax(S.costs) #find max cost and position
        if max_cost > threshold: #if max cost exceeds threshold declare change
            change_locations.append(dtime - w + position)
            detection_times.append(dtime)
    #above is the time of change and location in data stream
    
#One issue with Example 4 is that we return the same change multiple times. 
#Instead we can start from the next w new points post detection, as in Example 5:
    
np.random.seed(1)
X = np.concatenate([np.random.normal(0, 1, 50), np.random.normal(5,5,150),
                   np.random.normal(20, 1, 100), np.random.poisson(3, 100)])
X = np.ndarray.tolist(X)

w = 50 #window size
K = 3 #quantiles to use
S = NUNC.nunc_local_state(X[0], w) #initial state
change_locations = [] #for storing changepoints and detection times
detection_times = []

threshold = 12
dtime = 1 #number of points observed
restart_timer = 0
for x in X[1:]:
    S.update(K, x) #update state with new point
    dtime += 1
    restart_timer -= 1
    if len(S.window) == w and restart_timer <= 0:
        (max_cost, position) = argmax(S.costs) #find max cost and position
        if max_cost > threshold: #if max cost exceeds threshold declare change
            change_locations.append(dtime - w + position)
            detection_times.append(dtime)
            restart_timer = w - 1 #test from w points post detection
#NB: we could adapt this to instead start w points post change location, but
#this would often still give duplicates

