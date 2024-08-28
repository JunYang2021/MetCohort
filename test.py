import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import PchipInterpolator

input_x = np.array([4.51290846, 4.64421844, 4.77553368, 5.03815842, 5.16946983
                       , 5.30078554, 5.46969366, 6.12626028, 6.25757027, 6.38888645
                       , 9.8030405, 14.19114494, 15.11032867, 16.68608093, 19.31235886
                       , 22.51746178, 22.91140175, 25.45994949, 25.59127998, 27.29834175
                       , 27.95490837, 28.08621788, 36.28104019, 37.5941658, 40.92806625
                       , 41.05940628, 43.29170609, 44.07958221, 44.21090698, 44.34220886
                       , 53.92808151, 56.34318161, 58.05024338, 58.96943665, 64.74721527
                       , 68.60496521, 68.73627472, 68.86759186, 0.18371256, 134.57887268])
input_y = np.array([3.4629693, 3.98822141, 3.72559524, 3.65993774, 4.31649935
                       , 3.59428024, 3.72559524, 3.65993774, 3.72559524, 5.76094937
                       , 8.05893135, 10.29126072, 11.26301098, 12.97010422, 18.09131622
                       , 20.7688179, 21.29406929, 23.78901291, 19.66706848, 20.7688179
                       , 31.40517998, 21.29406929, 31.40517998, 33.76880646, 39.72927856
                       , 40.25452805, 41.88031387, 43.06212997, 43.19343948, 43.19343948
                       , 52.12042389, 54.4052887, 56.63759995, 56.50629425, 62.28408051
                       , 66.09214783, 66.48608398, 66.55174255, 0.18371256, 134.57887268])

lowess = sm.nonparametric.lowess(input_y, input_x, frac=1 / 5)

lowess_x = lowess[:, 0]
lowess_y = lowess[:, 1]

f_rt = interp1d(lowess_x, lowess_y, bounds_error=False, fill_value='extrapolate')
test_x = np.arange(0, 135, 0.1)
plt.plot(test_x, f_rt(test_x), color='green', linewidth=1)
# sorted_indices_rt = np.argsort(lowess_x)
# uni_x = lowess_x[sorted_indices_rt]
# uni_y = lowess_y[sorted_indices_rt]
# test_x = np.arange(0, 135, 0.1)
# pchip_rt = PchipInterpolator(uni_x, uni_y)
# plt.plot(test_x, pchip_rt(test_x), color='green', linewidth=1)

# Plot original data points
plt.scatter(input_x, input_y, label='Original Data', color='blue', alpha=0.5, s=1)

# Plot LOWESS trend
plt.scatter(lowess_x, lowess_y, label='LOWESS Trend', color='red', s=1)

# Add labels and title
plt.xlabel('Input X')
plt.ylabel('Input Y')
plt.title('LOWESS Smoothing')
plt.legend()

# Display the plot
plt.show()
