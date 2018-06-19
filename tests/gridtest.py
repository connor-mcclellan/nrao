from func import plot_grid, loadobj
import matplotlib.pyplot as plt

data_cube = loadobj('data_cube')
masks = loadobj('masks')
rejects = loadobj('rejects')
snr_vals = loadobj('snr_vals')
rlencatalog = loadobj('rangelencatalog')
    
def pick_handler(event):
    mouseevent = event.mouseevent
    artist = event.artist
    
def onpick(event):
    thisline = event.artist
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ind = event.ind
    points = tuple(zip(xdata[ind], ydata[ind]))
    print('onpick points:', points)

fig.canvas.mpl_connect('pick_event', onpick)

plot_grid(data_cube, masks, rejects, snr_vals, rlencatalog)
plt.suptitle('Interactivity Test')
plt.show()
