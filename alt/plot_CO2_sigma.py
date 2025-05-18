import numpy as np
import matplotlib.pyplot as plt

rho_c = (2440, 2240)    # density kg/m3 for min CO2 and max CO2
CO_c = (70, 90)         # min CO2 and max CO2
fcd_c = (16, 16)

rho_t = (470, 430)    # density kg/m3 for min CO2 and max CO2
CO_t = (137, 556)         # min CO2 and max CO2
fcd_t = (16, 16)

# Define the range of x values
x = np.linspace(1, 50, 100)

# Define the linear functions
y_upper_c = x/fcd_c[1]*CO_c[1]/rho_c[1]/1000
y_lower_c = x/fcd_c[0]*CO_c[0]/rho_c[0]/1000
y_upper_t = x/fcd_t[1]*CO_t[1]/rho_t[1]/1000
y_lower_t = x/fcd_t[0]*CO_t[0]/rho_t[0]/1000

# Plot the upper and lower limit functions
plt.plot(x/fcd_c[1], y_upper_c, color='red')
plt.plot(x/fcd_c[0], y_lower_c, color='red')

plt.plot(x/fcd_t[1], y_upper_t, color='blue')
plt.plot(x/fcd_t[0], y_lower_t, color='blue')

# Fill the area between the two functions
plt.fill_between(x/fcd_c[1], y_lower_c, y_upper_c, color='red', alpha=0.3)
plt.fill_between(x/fcd_t[1], y_lower_t, y_upper_t, color='blue', alpha=0.3)

# Add labels and title
plt.xlabel('x')
plt.ylabel('y')
plt.title('Area between Upper and Lower Limit Linear Functions')
plt.legend()

# Show the plot
plt.show()
