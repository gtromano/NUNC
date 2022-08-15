# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:56:29 2022

@author: austine
"""
from sortedcontainers import SortedList
from collections import deque
import math
    
class nunc_local_state:
    #this creates the initial sorted list and deque (window) for use with local
    #note we have length + 1 as we need to add the next point and start testing
    #when we use the nunc_local_update as per below
    
    '''
    Class for storing the window of data for use with nunc_local_update.
    
    Description
    -------------
    
    Instantiating this class takes inputted data and stores it in both an
    ordered tree and a deque. The deque is the window of points for use with
    NUNC Local, and the tree is that window of points in ordered search tree
    form so that quantiles can be computed in O(log w) time, where w is the 
    window size.
    
    The instantiation should takes the observations x, and stores them in a 
    window of size w. If the length of x is less than w (it only needs to be
     a single point) then the window is instantiated as a size w deque and 
    this can then be filled using the update() method.
    If the length of x exceeds w then we only store the most recent w points.
    
    The class comes with a built in method, update, that updates the window 
    and tree with a new observation and performs the NUNC Local algorithm on 
    this updated window. 
    
    The NUNC Local algorithm works by performing a search over the points 
    in the window for a change in distrbution. Each iteration of the NUNC 
    algorithm has complexity O(KW), where K is the number of quantiles and W
    is the window size.
    
    Input
    -------------
    x: list
        List of observations to store in the initial window. If the length
        exceeds w then only the last w observations are stored.
    w: int
        Positive integer for the size of the window.
    
    
    Attributes
    --------------
    tree: SortedList
        Ordered tree of observations in the window.
    window: deque
        Deque containing window of points. Initialises at length w.
    costs: List
        A list of costs computed by the NUNC algorithm. These costs correspond
        to the value of the test statistic for a change at each of the points
        in the window.
    
    Methods
    -------------
    
    update(x):
        Method to update the tree and window with a new observation, x. Once
        the window has been updated the NUNC algorithm is performed and then
        stored in the costs attribute.
    '''

    def __init__(self, x : list, w : int):
        
        try:
            list(x)
        except:
            try:
                if isinstance(x, float) or isinstance(x, int):
                    x = [x]
            except:    
                raise TypeError("x must be a listlike object")
        if not all(isinstance(y, (int, float)) for y in x):
            raise TypeError("x must be a list of floats / ints")
        if not (isinstance(w, int)):
            raise TypeError("w must be a positive integer")
        if w <= 0:
            raise ValueError("w must be a positive integer")   
        
        if len(x) > w: #only store at most the last w points
            x = x[-w:]
        
        self.tree = SortedList(x)
        self.window = deque(x, maxlen = w)
        self.costs = None
        
    def update(self, k : int, x : float) :
        
        '''
        Method to update the tree and window with a new observation. Once
        the window has been updated the NUNC algorithm is performed to update
        the costs. These costs correspond to the test statistic for a change
        at each point within the window.
        
        The function takes as input k, the number of quantiles to use for the
        NUNC algorithm; and x, the new point to update.
        
        Input
        ----------
        k: int
            Number of quantiles used by NUNC.
        x: float
            The next observation to store in the window. (Can also be type int 
            but will be converted.)

        '''
        
        if len(self.window) == self.window.maxlen :
            self.tree.remove(self.window[0])
        self.window.append(x)
        self.tree.add(x)
        #note we do not combine this logic as above because we fill the window
        #one iteration before we need to start removing points
        if len(self.window) == self.window.maxlen:
            self.costs = self.__nunc_local_update(k, x)
        
    def __nunc_local_update(self, k : int, x : float) : 
        
        '''
        Function for performing a single iteration of NUNC Local.
        
        Description
        ---------------
        
        This function performs a single step of the NUNC Local algorithm.
        That is, it tests a single window of data for a change in distribution.
        
        Once NUNC has been carried out on the window, a set of costs are 
        returnedalong with the updated window of data. In order to detect
        changes in the data stream, a user would check if the max of the
        returned costs exceed a given threshold. If so, then a change is 
        detected by NUNC. For more details see the examples.
                
        Input
        ----------

        k: int
            Number of quantiles used by NUNC.
        x: float
            The next observation to store in the window. (Can also be type int 
            but will be converted.)
        
        Returns
        ----------
        A list of costs of length w-1 denoting the test statistic at the points
        in the window. 
        '''
        
        if not (isinstance(k, int)):
            raise TypeError("K must be an integer")   
        if k <= 0:
            raise ValueError("K must be a positive integer")
        try:
            float(x)
        except:
            raise TypeError("x must be a float")
        
        def quantiles(data,k,w) :
        #this function computes k quantiles in data from a window of size w
        #the expression for the probabilities is taken from Zou (2014) and 
        #Haynes (2017) and is designed to add emphasis to the tails of the 
        #probability distribution. Both authors show doing so increases power.
        #data is the data to use for computing the quantiles
        #k is the number of quantiles and w is the window size
        

            def quantile(prob) :
                #this function works out the quantile value
                #it uses linear interpolation, ie quantile type = 7
                h = (len(data) - 1) * prob
                h_floor = int(h)
                if h_floor == h:
                    return data[h]
                else:
                    non_int_part = h - h_floor 
                    lower = data[h_floor]
                    upper = data[h_floor + 1]
                    return lower + non_int_part * (upper - lower)
            c = math.log(2*w-1)   #weight as in Zou (2014)
            probs = [(1/(1+(2*(w-1)*math.exp((-c/k)*(2*i-1)))))
                     for i in range(k)]
            return [quantile(p) for p in probs]    
        
        def eCDF_vals(data,quantile):    
            #used to return value of eCDF, not cost, at a given quantile  
            #data is the tree of data used to compute the ecdf
            #quantile is the numeric quantile value
            
            left = data.bisect_left(quantile)
            right = data.bisect_right(quantile)
            #value is number of points to left of quantile, plus 0.5 times
            #the points equal to the quantile
            val = (left+0.5*(right-left))/len(data)
            return val
        
        def one_point_emp_dist(data, quantile):
            #function for computing empirical CDF for data at a set quantile
            #ie, is a point less than, equal to, or greater than the quantile
            #data is an array of numerics
            #quantile is the quantile to evaluate the eCDF at.
            if(data < quantile):
                return(1)
            elif (data == quantile):
                return(0.5)
            else: 
                return(0)
            
        def cdf_cost(cdf_val, seg_len):
            #function for computing the likelihood function
            #cdf_val is the value of the eCDF at a set quantile
            #seg_len is the length of the data used 
            if(cdf_val <= 0 or cdf_val >= 1):
                return(0) #avoids rounding error, does not affect result
            conj = 1 - cdf_val
            cost = seg_len * (cdf_val * math.log(cdf_val)
                              - conj * math.log(conj))
            return(cost)
        
        def update_window_ecdf_removal(data_to_remove, quantiles, current_ecdf,
                                       current_len):
            #this function takes a set of K current eCDFs values
            #computed from a data set of length current_len, and removes a 
            #point from this data. 
            #The function then returns the updated CDF values following removal
            #of this point.
            num_quantiles = len(quantiles)
            for i in range(num_quantiles):
                current_ecdf[i] *= current_len
                current_ecdf[i] -= one_point_emp_dist(data_to_remove,
                                                      quantiles[i]) 
                current_ecdf[i] /= (current_len - 1)  
            return current_ecdf
        
        tree = self.tree #extract tree and window of points
        window = self.window
        w = len(window)  #compute window size from inputted data
        Q = quantiles(tree,k,w) #update quantiles
        full_cdf_vals = [eCDF_vals(tree, q) for q in Q] #full data eCDF
        right_cdf_vals = full_cdf_vals.copy() #updates as we search for change 
        full_cost = sum(cdf_cost(val, w) for val in full_cdf_vals) 
        segment_costs = list() #used for storing costs of segmented window
        current_len = w #current length of right segment, updates iteratively
        left_cdf_vals = [0] * len(Q) #update as we search window for points
        for i in range(0, w-1): #window updates are O(K)
        #as we loop over the window we "move" points from the right to
        #left segment, and update the eCDFs. This provides an O(K) cost
        #for updating the eCDFs for each segment.
            right_cdf_vals = update_window_ecdf_removal(window[i], Q,
                                                right_cdf_vals, current_len)
            #remove points from RHS iteratively and update eCDF
            current_len -= 1
            for j in range(len(Q)): #update LHS using RHS and full eCDFs
                left_cdf_vals[j] = ((full_cdf_vals[j]*w - right_cdf_vals[j]*current_len) / (w - current_len))
            #compute costs of segmented data
            left_cost = sum([cdf_cost(val, w - current_len)
                             for val in left_cdf_vals])
            right_cost = sum([cdf_cost(val, current_len) 
                              for val in right_cdf_vals])
            segment_costs.append(left_cost + right_cost)
        #update cost function
        costs = [2*(cost - full_cost) for cost in segment_costs]
        return costs

