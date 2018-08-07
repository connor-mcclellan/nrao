from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
from collections import OrderedDict

def colorcolor_plot(xfreq, x_flux, x_err, yfreq, y_flux, y_err, numfreq, num_flux, num_err):
    xgood = [x_flux[i]>3.*x_err[i] for i in range(len(x_flux))]
    ygood = [y_flux[i]>3.*y_err[i] for i in range(len(y_flux))]
    good = [a and b for a, b in zip(xgood, ygood)]
    plotx = num_flux[good]/x_flux[good]
    ploty = num_flux[good]/y_flux[good]
    plotx_err = x_flux[good]*np.sqrt((x_err[good]/x_flux[good])**2 + (num_err[good]/num_flux[good])**2) #error is fractional error added in quad
    ploty_err = x_flux[good]*np.sqrt((y_err[good]/y_flux[good])**2 + (num_err[good]/num_flux[good])**2)
    alphas = [0., 1., 2., 3.]
    for i, alpha in enumerate(alphas):
        plt.vlines(((float(numfreq.split('GHz')[0])/float(xfreq.split('GHz')[0]))**alpha), 0, 100, linestyles='dashed', colors=colors[i], label='alpha={}'.format(alpha))
        plt.hlines(((float(numfreq.split('GHz')[0])/float(yfreq.split('GHz')[0]))**alpha), 0, 100, linestyles='dashed', colors=colors[i], label='alpha={}'.format(alpha))
    plt.errorbar(plotx, ploty, xerr=plotx_err, yerr=ploty_err, linestyle='', marker='.', color='k', label='Fluxes')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Log {num}/{x} Flux'.format(num=numfreq,x=xfreq))
    plt.ylabel('Log {num}/{y} Flux'.format(num=numfreq,y=yfreq))
    plt.xlim([1e-2, 1e2])
    plt.ylim([5e-1, 1e2])
    handles, labels = plt.gca().get_legend_handles_labels()
    label = OrderedDict(zip(labels, handles))
    plt.legend(label.values(), label.keys())
    plt.title('Color v. Color for Sources in W51 e2')
    #plt.show()
    plt.savefig('/users/bmcclell/nrao/documentation/colorcolor/colorcolor_e1e2_{}div{}{}.png'.format(numfreq, xfreq, yfreq), overwrite=True)
    

table = Table.read('/users/bmcclell/nrao/cat/45-93-226GHz_photometered_adjustedRADEC.dat', format='ascii')
ind = [table['rejected']==0]
xflux = table['226.1GHz_Ellipse_sum'][ind]
xerr = table['226.1GHz_Annulus_rms'][ind]
yflux = table['45.0GHz_Ellipse_sum'][ind]
yerr = table['45.0GHz_Annulus_rms'][ind]
numflux = table['93.0GHz_Ellipse_sum'][ind]
numerr = table['93.0GHz_Annulus_rms'][ind]
colorcolor_plot('226.1GHz', xflux, xerr, '45.0GHz', yflux, yerr, '93.0GHz', numflux, numerr)


