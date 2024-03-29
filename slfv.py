# We simulate the spatial Lambda-Fleming-Viot model

# This line is useless.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import colors
from IPython import embed
from os import path
import sys

class poissonpp:
	"""
	We sample a Poisson point process
	"""
	def __init__(self, lam, space_horizon, time_horizon):
		# Initial values
		self.lam = lam
		self.space_horizon = space_horizon
		self.time_horizon = time_horizon

	def sample(self):

		# Creates a point_num x 2 array, with first column ordered times and second column associated space position

		self.point_num = np.random.poisson(self.lam*self.space_horizon*self.time_horizon,size = None)
		self.points = np.random.rand(self.point_num, 2)
		self.points = self.points*[self.time_horizon, self.space_horizon]
		self.points = my_sort(self.points)
		
		return self.points

def my_sort(my_matrix):

	# Automatically sorts along the first column (in our case according to times)
	indexes = np.argsort(my_matrix, axis = 0)[:,0]

	return my_matrix[indexes, :]


class slfv:

	"""
	We simulate the evolution of a spatial Lambda-Fleming-Viot process
	"""

	def __init__(self, slfv_init, space_points, dx, space_horizon):

		# Initial state of the system
		self.state = slfv_init
		self.space_points_len = len(self.state)
		self.space_points = space_points
		self.dx = dx

		# Other parameters (mainly for PPP)
		self.space_horizon = space_horizon
		self.time_horizon = 10.0
		self.lam = 0.07

		# Counters
		self.count = 0
		self.cur_time =0

		# Parameters of the SLFV
		self.impact = 0.1
		self.radius = 3.0

		if self.radius < self.dx:
			print( "\n\n ERROR: The radius of the impact zones is smaller than the spatial scale \n\n")


		# We sample the poisson times of the system 
		self.ppp =poissonpp(self.lam, self.space_horizon, self.time_horizon)
		self.jump_points=self.ppp.sample()
		self.available_points = self.ppp.point_num

		self.next_jump = self.jump_points[0,:]

		# Variables that we need for the jumps
		self.local_average = 0
		self.local_average_count =0
		self.choice =0

	def go_to_time(self, next_time):

		self.next_time = next_time

		while self.next_jump[0] < self.next_time:
			while (self.next_jump[0] < self.next_time and self.count<self.available_points-1):

				self.jump(self.next_jump[1])

				# We adjourn all variables
				self.count += 1
				self.cur_time = self.next_jump[0]
				self.next_jump = self.jump_points[self.count, :]


			# If we exit because we don't have enough points we sample new points
			if self.count >= self.available_points-1:				

				self.jump_points=self.ppp.sample() + [self.cur_time,0]
				self.available_points = self.ppp.point_num
				self.count = 0
					


	def jump(self, loc):
		
		# First we choose the direction of the reproduction event
		self.average(loc)
		self.choice = np.random.binomial(1,self.local_average)

		# Growth
		if self.choice == 1:
		
			for i in range(self.space_points_len):
				if np.abs(self.space_points[i]-loc) < self.radius:
					
					self.state[i] = self.state[i]+self.impact*(1-self.state[i])

		# Decrease
		if self.choice ==0:

			for i in range(self.space_points_len):
				if np.abs(self.space_points[i]-loc) < self.radius:
					
					self.state[i] = self.state[i]*(1-self.impact)




	def average(self, loc):

		# This function builds the local average of the process around "loc"
		self.local_average=0
		self.local_average_count=0.0
		for i in range(self.space_points_len):
			if np.abs(self.space_points[i]-loc) < self.radius:
				self.local_average += self.state[i]
				self.local_average_count+=1.0
		
		# We normalize the result
		self.local_average = self.local_average/self.local_average_count


# This does the animation

def animate(i):
	# Real time is:

	ani_time = delta_t*(i+0.07*(i**2))

	# Redefine the plot
	my_lines.set_data(space_points, my_slfv.state)

	# Set the new time
	time_text.set_text("Time = {:2.3f}".format(ani_time) )
	# We print the step we are in:
	sys.stdout.flush()
	#sys.stdout.write("\r Step = {}, Value = {}, Derivative = {}".format(i, fkpp_sample.state_a[middle], (delta_x**(-2))*(fkpp_sample.state_a[middle-1]+fkpp_sample.state_a[middle+1]-2*fkpp_sample.state_a[middle])))
	sys.stdout.write("\r Step = {}".format(i))
	# And we do the next step:
	my_slfv.go_to_time(ani_time)
	return [my_lines,] + [time_text,]

# We define the parameters of the SLFV

# Parameters of the domain
space_horizon = 90.0
dx = 0.05
space_points = np.arange(0,space_horizon,dx)
space_len = len(space_points)

# Initial condition
slfv_init = 0.5 + np.zeros(shape = space_len)

# We initialize the SLFV process
my_slfv = slfv(slfv_init, space_points, dx, space_horizon)

#We set up the picture
fig       = plt.figure()
ax        = plt.axes(xlim=(0, space_horizon-5.5), ylim = (-0.02, 1.02))
time_text = ax.text(0.05, 0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)
my_lines,  = ax.plot([],[], lw = 2)
plt.title("The neutral SLFV process")

def init():
    my_lines.set_data([], [])
    return my_lines,

# We let the animation go.
delta_t=0.005
ani = FuncAnimation(fig, animate, init_func=init, frames=2000, interval=10, blit=True)
ani.save('slfv.gif', writer='imagemagick')
# ani.save(filename = 'neutral_slfv.html')

# ani       = animation.FuncAnimation(fig, animate, frames=400, interval = 70, blit = True)

# ani.save(filename = 'neutral_slfv.html', extra_args=['-vcodec', 'libx264'], bitrate = 20000)


# INSTRUCTION FOR PUTTING VIDEO IN PRESENTATION.

# 1) RUN: ffmpeg -i <input> -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 output.mp4
#	 on powershell. The result (output.mp4) is the video you will use.
# 2)  IN Latex, with package movie9 write:
#   \includemedia[
#  width=0.7\linewidth,
#  totalheight=0.7\linewidth,
#  activate=onclick,
#  %passcontext,  %show VPlayer's right-click menu
#  addresource=ballistic_out.mp4,
#  flashvars={
#    %important: same path as in `addresource'
#    source=ballistic_out.mp4
#  }
#]{\fbox{Click!}}{VPlayer.swf}


# # Definition of the Poisson PP
# ppp = poissonpp(lam, space_horizon,time_horizon)

# # To check that the Poisson point process works
# points = ppp.sample()
# plt.scatter(points[:,0], points[:,1])
# plt.show()








