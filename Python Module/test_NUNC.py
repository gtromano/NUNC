# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 13:52:52 2022

@author: austine
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 09:42:45 2022

@author: austine
"""
import NUNC
import pytest
import numpy as np

def test_nunc_local_state1():
    np.random.seed(1)
    X = np.random.normal(0, 1, 20)
    S = NUNC.nunc_local_state(X[:19], len(X))
    K = 1
    
    assert (S.window[0] == X[0])
    assert(len(S.window) == len(X)-1)
    assert(S.tree[0] == min(X[:19]))
    assert(S.tree[18] == max(X[:19]))
    assert(S.costs is None)
    
    S.update(K, X[19])
    assert(S.window[0] == X[0])
    assert(S.window[19] == X[19])
    assert(len(S.window) == len(X))
    assert(max(S.costs) == 2.2820159898644965)

def test_nunc_local_state2():
    np.random.seed(2)
    X = np.random.normal(-2, 4, 40)
    S = NUNC.nunc_local_state(X[:39], len(X))
    K = 4
    
    assert (S.window[0] == X[0])
    assert(len(S.window) == len(X)-1)
    assert(S.tree[0] == min(X[:39]))
    assert(S.tree[38] == max(X[:39]))
    assert(S.costs is None)
    
    S.update(K, X[39])
    assert(S.window[0] == X[0])
    assert(S.window[39] == X[39])
    assert(len(S.window) == len(X))
    assert(max(S.costs) == 4.133634955432745)

def test_nunc_local_state_update1():
    np.random.seed(1)
    X = np.random.normal(-2, 4, 40)
    S = NUNC.nunc_local_state(X[:40], len(X) + 1)
    y = [1.2, 2.1]
    K = 2
    
    S.update(K, y[0])
    assert(max(S.costs) == 3.7805882252859977)
    S.update(K, y[1])
    assert(max(S.costs) == 4.222114269620466)
    
    assert(S.window[0] == X[1])
    assert(S.window[39] == y[0])
    assert(S.window[40] == y[1])
    assert(X[0] not in S.tree)
    
def test_nunc_local_state_update2():
    np.random.seed(2)
    X = np.random.normal(0, 1, 40)
    S = NUNC.nunc_local_state(X[:40], len(X) + 1)
    y = [-9, 0]
    K = 5
    
    S.update(K, y[0])
    assert(S.costs[0] == -0.4684940490981351)
    S.update(K, y[1])
    assert(min(S.costs) == 0.4826754106270936)
    
    assert(S.window[0] == X[1])
    assert(S.window[39] == y[0])
    assert(S.window[40] == y[1])
    assert(X[0] not in S.tree)

def test_nunc_local_state_update3():
    np.random.seed(1)
    X = np.random.normal(0, 1, 10)
    S = NUNC.nunc_local_state(X[:9], len(X))
    K = 1
    
    S.update(K, X[9])
    assert(S.costs[5] == 0.9483775336306544)
    
    Y = np.random.poisson(1, 100)
    for y in Y:
        S.update(K, y)
        assert(S.window[9] == y)
    assert(S.costs[0] == 0.4186175690251659)
    
def test_nunc_local_state_update4():
    np.random.seed(1)
    X = np.random.normal(0, 1, 10)
    S = NUNC.nunc_local_state(X, 5)
    
    assert(len(S.window) == 5)
    assert(S.window[0] == X[5])
    assert(S.window[4] == X[9])
    
def test_nunc_local_update5():
    np.random.seed(1)
    X = np.random.normal(0, 1, 1)
    S = NUNC.nunc_local_state(X, 5)
    
    assert(S.window[0] == X[0])
    assert(S.window.maxlen == 5)
    
def test_nunc_local_state_update6():
    np.random.seed(1)
    X = np.random.normal(0, 1, 10)
    S = NUNC.nunc_local_state(X[0], len(X))
    K = 1
    
    for x in X[1:]:        
        S.update(K, x)
    assert(S.costs[5] == 0.9483775336306544)
    
    Y = np.random.poisson(1, 100)
    for y in Y:
        S.update(K, y)
        assert(S.window[9] == y)
    assert(S.costs[0] == 0.4186175690251659)
    
def test_nunc_local_state_x_TypeError1():
    with pytest.raises(TypeError):
        NUNC.nunc_local_state("Type Error", 5)
    
def test_nunc_local_state_x_TypeError2():
    with pytest.raises(TypeError):
        NUNC.nunc_local_state(["Type Error"], 5)
    
def test_nunc_local_state_x_TypeError3():
    with pytest.raises(TypeError):
        NUNC.nunc_local_state([1, "Type Error"], 5)

def test_nunc_local_state_w_TypeError1():
    np.random.seed(1)
    X = np.random.normal(0, 1, 1)
    with pytest.raises(TypeError):
        NUNC.nunc_local_state(X, "TypeError")
        
def test_nunc_local_state_w_TypeError2():
    np.random.seed(1)
    X = np.random.normal(0, 1, 1)
    with pytest.raises(TypeError):
        NUNC.nunc_local_state(X, 3.0)
        
def test_nunc_local_state_w_ValueError1():
    np.random.seed(1)
    X = np.random.normal(0, 1, 1)
    with pytest.raises(ValueError):
        NUNC.nunc_local_state(X, 0)
        
def test_nunc_local_state_w_ValueError2():
    np.random.seed(1)
    X = np.random.normal(0, 1, 1)
    with pytest.raises(ValueError):
        NUNC.nunc_local_state(X, -1)

def test_nunc_local_state_k_TypeError():
    X = np.random.normal(0, 1, 10)
    S = NUNC.nunc_local_state(X, len(X))
    
    with pytest.raises(TypeError):
        S.update("TypeError", 1.2)
        
def test_nunc_local_state_k_ValueError():
    X = np.random.normal(0, 1, 10)
    S = NUNC.nunc_local_state(X, len(X))
    
    with pytest.raises(ValueError):
        S.update(-2, 1.2)
    
def test_nunc_local_state_x_TypeError():
    X = np.random.normal(0, 1, 10)
    S = NUNC.nunc_local_state(X, len(X))
    
    with pytest.raises(TypeError):
        S.update(2, "TypeError")

pytest.main()
