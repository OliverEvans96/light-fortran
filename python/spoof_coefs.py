# File Name: spoof_coefs.py
# Description: Make up absorption/scattering coefficients for RTE2D
# Created: Mon Jan 16, 2017 | 04:52pm EST
# Last Modified: Mon Jan 16, 2017 | 05:22pm EST

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#                           GNU GPL LICENSE                            #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#                                                                      #
# Copyright Oliver Evans 2017 <oliverevans96@gmail.com>                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
#                                                                      #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# Imports
from numpy import *

# Box dimensions
imax = 50
jmax = 50
xmin = 0
xmax = 10
ymin = 0
ymax = 10

# Allocate arrays (exclude upper endpoint)
xx = linspace(xmin,xmax,imax+1)[:-1]
yy = linspace(ymin,ymax,imax+1)[:-1]
xx,yy = meshgrid(xx,yy,indexing='ij')

# Function
def prob_dist(xx,yy):
    return exp(-((xx-xmin)/(xmax-xmin))**2+(yy-ymin)/(ymax-ymin))

# Write array
arr = prob_dist(xx,yy)
savetxt('../data/coefs/spoof_prob_dist2d.txt',arr,fmt='%10.3E',delimiter='')
