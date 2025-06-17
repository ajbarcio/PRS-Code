from modules.utils import line_of_action
import numpy as np

center, force = line_of_action([1,1,1], [[-np.sqrt(2)/2, np.sqrt(2)/2]])
print(center, force)
