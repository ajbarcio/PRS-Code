from modules.utils import line_of_action
import numpy as np

def test_line_of_action():
    center, force = line_of_action([np.sqrt(2)/2,np.sqrt(2)/2,1], [-np.sqrt(2)/2, np.sqrt(2)/2])
    print(center, force)
