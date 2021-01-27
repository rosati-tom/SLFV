# We simulate the 2D spatial Lambda-Fleming-Viot model

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

		# Creates a point_num x 3 array, with first column ordered times and second column associated space position

		self.point_num = np.random.poisson(self.lam*(self.space_horizon**2)*self.time_horizon,size = None)
		self.points = np.random.rand(self.point_num, 3)
		self.points = self.points*[self.time_horizon, self.space_horizon, self.space_horizon]
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
		self.time_horizon = 20.0
		self.lam = 0.1

		# Counters
		self.count = 0
		self.cur_time =0

		# Parameters of the SLFV
		self.impact = 0.2
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

		# This gives the maximum  distance the index has to run through:
		self.int_limit = int(self.radius/dx)+10

	def go_to_time(self, next_time):

		self.next_time = next_time

		while self.next_jump[0] < self.next_time:
			while (self.next_jump[0] < self.next_time and self.count<self.available_points-1):
				
				self.next_jump_loc = self.next_jump[1:3]
				self.jump()

				# We adjourn all variables
				self.count += 1
				self.cur_time = self.next_jump[0]
				self.next_jump = self.jump_points[self.count, :]


			# If we exit because we don't have enough points we sample new points
			if self.count >= self.available_points-1:				

				self.jump_points=self.ppp.sample() + [self.cur_time,0]
				self.available_points = self.ppp.point_num
				self.count = 0
					


	def jump(self):

		# We associate two indices to the locations:
		self.int_loc = self.next_jump_loc/dx

		self.int_loc_bot = np.where(self.int_loc - self.int_limit>0, self.int_loc - self.int_limit, 0)
		self.int_loc_bot = [int(self.int_loc_bot[0]), int(self.int_loc_bot[1])] 

		self.int_loc_top = np.where(self.int_loc + self.int_limit < self.space_points_len, self.int_loc + self.int_limit, self.space_points_len)
		self.int_loc_top = [int(self.int_loc_top[0]), int(self.int_loc_top[1])] 

		# print( self.int_limit, self.int_loc, self.next_jump_loc, self.int_loc_top, self.int_loc_bot)

		
		# First we choose the direction of the reproduction event
		self.average()
		self.choice = np.random.binomial(1,self.local_average)

		# Growth
		if self.choice == 1:
		
			for i in range(self.int_loc_bot[0], self.int_loc_top[0]):
				for j in range(self.int_loc_bot[1], self.int_loc_top[1]):
						
					self.dist = np.linalg.norm(self.space_points[[i,j]] - self.next_jump_loc)
					
					if self.dist < self.radius:

						self.state[i,j] = self.state[i,j]+self.impact*(1-self.state[i,j])

		# Decrease
		if self.choice ==0:
			
			for i in range(self.int_loc_bot[0], self.int_loc_top[0]):
				for j in range(self.int_loc_bot[1], self.int_loc_top[1]):
						
					self.dist = np.linalg.norm(self.space_points[[i,j]] - self.next_jump_loc)
					
					if self.dist < self.radius:

						self.state[i,j] = self.state[i,j]*(1-self.impact)




	def average(self):

		# This function builds the local average of the process around "loc"
		self.local_average=0
		self.local_average_count=0.0

		for i in range(self.int_loc_bot[0], self.int_loc_top[0]):
			for j in range(self.int_loc_bot[1], self.int_loc_top[1]):
					
				self.dist = np.linalg.norm(self.space_points[[i,j]] - self.next_jump_loc)
				
				if self.dist < self.radius:

					self.local_average += self.state[i,j]
					self.local_average_count+=1.0


		# We normalize the average
		self.local_average = self.local_average/self.local_average_count


# This does the animation
def animate(i):
	# Real time is:

	ani_time = delta_t*(i+0.05*(i**2))
	my_slfv.radius =np.maximum(0.2, my_slfv.radius-0.01*ani_time)

	# Redefine the plot
	my_im.set_data(my_slfv.state)


	# Set the new time
	time_text.set_text("Time = {:2.3f}".format(ani_time) )
	# We print the step we are in:
	sys.stdout.flush()
	#sys.stdout.write("\r Step = {}, Value = {}, Derivative = {}".format(i, fkpp_sample.state_a[middle], (delta_x**(-2))*(fkpp_sample.state_a[middle-1]+fkpp_sample.state_a[middle+1]-2*fkpp_sample.state_a[middle])))
	sys.stdout.write("\r Step = {}".format(i))
	# And we do the next step:
	my_slfv.go_to_time(ani_time)
	return [my_im,] + [time_text,]

# We define the parameters of the SLFV

# Parameters of the domain
space_horizon = 60.0
dx = 0.04
space_points = np.arange(0,space_horizon,dx)
space_len = len(space_points)

# Initial condition
slfv_init = 0.5 + np.zeros(shape = (space_len, space_len))

# We initialize the SLFV process
my_slfv = slfv(slfv_init, space_points, dx, space_horizon)

#We set up the picture
fig       = plt.figure()
ax        = plt.axes()
plt.axis('off')
# ax        = plt.axes(xlim=(0, space_horizon-5.5), ylim = (0, space_horizon-5.5))
time_text = ax.text(0.05, 0.95,'',horizontalalignment='left',verticalalignment='top', transform=ax.transAxes)

my_im     = ax.imshow(slfv_init, interpolation='none', vmin = 0, vmax = 1, cmap ='bwr')
# Change color map with cmap setting

plt.title("The neutral SLFV process in two dimensions")

# Picture for the particle system
# im1       = ax1.imshow(brwre_histo*(B_NUM_EXTRA/DIM_BOX)**2, interpolation='none', origin='low', vmin = 0, vmax = 0.034*intensity*B_NUM_EXTRA,  \
# 	extent=[space_histo_x[0]/DIM_BOX, space_histo_x[-1]/DIM_BOX, space_histo_y[0]/DIM_BOX, space_histo_y[-1]/DIM_BOX], cmap = plt.get_cmap('jet'))

embed()

def init():

    my_im.set_data(my_slfv.state)
    return my_im,

# We let the animation go.
delta_t=0.005
ani = FuncAnimation(fig, animate, init_func=init, frames=500, interval=2, blit=True)
ani.save('2D_slfv.gif', writer='imagemagick')
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








