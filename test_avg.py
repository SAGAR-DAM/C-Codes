# import AvgModule
import numpy as np
from CPP_Module import AvgModule

# Example NumPy array
arr = np.random.uniform(low=0, high=30,size=300)

# Call the average function from the C++ module
result = AvgModule.cpp_average(arr)

print(f"The average of the array is: {result}")

# print(dir(AvgModule))