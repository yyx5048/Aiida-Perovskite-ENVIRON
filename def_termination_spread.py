def find_termination_atomicpsread(c,spread):
    """
    read the calculation label and return the termiantion atoms
    and for atomispread index, which return the atomicspread list
    """
    all_atm=c.inp.structure.get_kind_names() # elements in alphametic order
    term_atm = c.label.split()[4] # element at the termination

    oxy_index = all_atm.index("O")
    term_index = all_atm.index(term_atm)

    default_spread = [0.5,0.5,0.5]
    default_spread[oxy_index] = spread
    default_spread[term_index] = spread

    return default_spread # return the atomicspread that only push away the termiantion cavaty.

