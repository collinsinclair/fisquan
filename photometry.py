import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from photutils import CircularAperture, aperture_photometry
from photutils.aperture import CircularAnnulus, CircularAperture
from photutils.background import SExtractorBackground
from photutils.detection import DAOStarFinder
from photutils.utils import calc_total_error


def find_sources(filepath, plot=False):
    header = fits.open(filepath)[0].header  # get header
    data = fits.open(filepath)[0].data  # get data
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)  # get stats

    daofind = DAOStarFinder(fwhm=7.35, threshold=10*std)  # DAOStarFinder object
    sources = daofind(data - median)  # find sources
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  
    print(sources)

    if plot:
        plt.close()
        fig, ax = plt.subplots(1, 1)
        img = ax.imshow(data, 'gist_gray', origin='lower',
                        vmin=mean-(3*std), vmax=mean+(3*std))  # plot image

        aps = CircularAperture(positions=zip(
            sources["xcentroid"], sources["ycentroid"]), r=7)
        aps.to_sky(wcs=WCS(header))  # convert to sky coordinates
        aps.plot(color='red')  # plot apertures

        cbr = plt.colorbar(img, ax=ax)
        cbr.set_label('cts/pix')
        ax.set_title(filepath[8:13])
        plt.show()

    return sources


def main():
    file_paths = ['masters/m34_B_master.fits',
                  'masters/m34_V_master.fits',
                  'masters/m37_B_master.fits',
                  'masters/m37_V_master.fits']
    sources = []
    for path in file_paths:
        sources.append(find_sources(path, True))


if __name__ == '__main__':
    main()
