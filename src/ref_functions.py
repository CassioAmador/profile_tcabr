"""
Function: Common functions and constants associated with Reflectometry.
TODO: ?
"""

k=8.9786629998410721
#(1/k)^2
k_2=0.012404425565580237

def freq2den(freq):
    """if frequency in Hz, returns den in m^-3"""
    return freq*freq*k_2

def den2freq(den):
    """if den in m^-3, returns frequency in Hz"""
    return k*p.sqrt(den)
