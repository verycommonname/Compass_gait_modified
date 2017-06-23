import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def update_line(num, data, line):
    line.set_ydata([0, data[num]])
    print ('line', line)
    print ('type', type(line))
    print ('line,', line,)
    return line, 

fig1 = plt.figure()
leg1, = plt.plot([0, 1], [0.2, 0.5], 'r-')

plt.xlim(0, 1)
plt.ylim(0, 2)
plt.xlabel('x')
plt.title('leg test')

data = np.linspace(0, 1.5)

line_ani = animation.FuncAnimation(fig1, update_line, 50, fargs= (data, leg1), 
                                   interval = 50, blit = True)
# line_ani.save('legs_video.mp4')

plt.show()
