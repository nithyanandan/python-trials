import numpy as NP
from astropy.io import fits
from astropy.io import ascii
import scipy.constants as FCNST
# matplotlib.use('TkAgg')
import matplotlib.pyplot as PLT
import matplotlib.animation as MOV
import geometry as GEOM
import interferometry as RI
import catalog as CTLG

catalog_file = '/Users/t_nithyanandan/Downloads/mwacs_b1_131016.csv'

catdata = ascii.read(catalog_file, data_start=1, delimiter=',')

dec_deg = catdata['DEJ2000']
ra_deg = catdata['RAJ2000']
fpeak = catdata['S150_fit']
ferr = catdata['e_S150_fit']
spindex = catdata['Sp+Index']
freq_catalog = 0.150 # in GHz
freq_resolution = 40.0 # in kHz
nchan = 256
chans = freq_catalog + (NP.arange(nchan) - 0.5 * nchan) * freq_resolution * 1e3 / 1e9
bpass = 1.0*NP.ones(nchan)

ctlgobj = CTLG.Catalog(freq_catalog, NP.hstack((ra_deg.reshape(-1,1), dec_deg.reshape(-1,1))), fpeak)

skymod = CTLG.SkyModel(ctlgobj)

A_eff = 16.0 * (0.5 * FCNST.c / (freq_catalog * 1e9))**2

intrfrmtr = RI.Interferometer('B1', [1000.0, 0.0, 0.0], chans, telescope='mwa',
                              latitude=-26.701, A_eff=A_eff, freq_scale='GHz')

Tsys = 440.0 # in Kelvin
t_snap = 20 * 60.0 # in seconds
# ha_range = 15.0*NP.arange(-1.0, t_snap/3.6e3, 1.0)
n_snaps = 16
lst_obs = (0.0 + (t_snap / 3.6e3) * NP.arange(n_snaps)) * 15.0 # in degrees
ha_obs = 1.0 * NP.zeros(n_snaps)
dec_obs = intrfrmtr.latitude + NP.zeros(n_snaps)

for i in xrange(n_snaps):
    intrfrmtr.observe(str(lst_obs[i]), Tsys, bpass, [ha_obs[i], dec_obs[i]], skymod, t_snap, fov_radius=30.0, lst=lst_obs[i])

intrfrmtr.delay_transform()

# fig1 = PLT.figure(figsize=(14,8))
# ax1 = PLT.axes([0.12,0.1,0.8,0.8])
# PLT.xlim(1e6*NP.min(intrfrmtr.lags)-1.0, 1e6*NP.max(intrfrmtr.lags)+1.0)
# PLT.ylim(0.5*NP.min(NP.abs(intrfrmtr.skyvis_lag)),2.0*NP.max(NP.abs(intrfrmtr.skyvis_lag)))
# PLT.xlabel(r'$\eta$ [$\mu$s]', fontsize=18)
# PLT.ylabel('Amplitude [Jy Hz]', fontsize=18)
# PLT.title('Delay Spectrum', fontsize=18)

# phase_center = ax1.text(0.05, 0.9, '', transform=ax1.transAxes, fontsize=18)

# l, = PLT.semilogy([], [], 'k+', markersize=10)

# def init():
#     l.set_xdata([])
#     l.set_ydata([])
#     phase_center.set_text('')
#     return l, phase_center

# def update(i, interferometer, line, pc_text):
#     line.set_xdata(1e6 * interferometer.lags)
#     line.set_ydata(NP.abs(interferometer.skyvis_lag[i,:]))
#     line.set_marker('+')
#     line.set_color('k')
#     line.set_markersize(10)
#     label_str = 'RA = {0:+.3f} deg, DEC = {1:+.2f} deg'.format(float(interferometer.timestamp[i])-interferometer.pointing_center[i,0], interferometer.pointing_center[i,1])
#     pc_text.set_text(label_str)
#     return line, pc_text

# anim = MOV.FuncAnimation(fig1, update, fargs=(intrfrmtr, l, phase_center), frames=intrfrmtr.vis_lag.shape[0], interval=400, blit=True, init_func=init)
# PLT.show()


# fig2 = PLT.figure(figsize=(14,8))
# ax2 = PLT.axes([0.12,0.1,0.8,0.8])
# PLT.xlim(NP.min(skymod.catalog.location[:,0]), NP.max(skymod.catalog.location[:,0]))
# PLT.ylim(NP.min(skymod.catalog.location[:,1])-5.0, NP.max(skymod.catalog.location[:,1])+5.0)
# PLT.xlabel(r'$\alpha$ [degrees]', fontsize=18)
# PLT.ylabel(r'$\delta$ [degrees]', fontsize=18)
# PLT.title('Sky Model', fontsize=18)
# PLT.grid(True)

# # l, = PLT.plot([], [], 'k.', markersize=10)
# phase_center = ax1.text(0.05, 0.9, '', transform=ax1.transAxes, fontsize=18)
# l, = PLT.plot(skymod.catalog.location[:,0], skymod.catalog.location[:,1], 'k.', markersize=1)

# def init():
#     l.set_xdata(skymod.catalog.location[:,0])
#     l.set_ydata(skymod.catalog.location[:,1])
#     l.set_marker('.')
#     phase_center.set_text('')
#     return l, phase_center

# def update(i, interferometer, line, pc_text):
#     # ra_center = float(interferometer.timestamp[i])-interferometer.pointing_center[i,0]
#     # dec_center = interferometer.pointing_center[i,1]
#     line.set_xdata(skymod.catalog.location[NP.asarray(interferometer.obs_catalog_indices[i]),0])
#     line.set_ydata(skymod.catalog.location[NP.asarray(interferometer.obs_catalog_indices[i]),1])
#     line.set_marker('+')
#     line.set_markersize(3)
#     line.set_color('r')
#     label_str = 'RA = {0:+.3f} deg, DEC = {1:+.2f} deg'.format(float(interferometer.timestamp[i])-interferometer.pointing_center[i,0], interferometer.pointing_center[i,1])
#     pc_text.set_text(label_str)
#     return line, pc_text

# anim = MOV.FuncAnimation(fig2, update, fargs=(intrfrmtr, l, phase_center), frames=intrfrmtr.vis_lag.shape[0], interval=400, blit=True, init_func=init)
# PLT.show()

fig = PLT.figure(figsize=(14,14))
ax1 = fig.add_subplot(211)
# fig, (ax1, ax2) = PLT.subplots(2,1,figsize=(14,12))
ax1.set_xlabel(r'$\eta$ [$\mu$s]', fontsize=18)
ax1.set_ylabel('Amplitude [Jy Hz]', fontsize=18)
ax1.set_title('Delay Spectrum', fontsize=18, weight='semibold')
ax1.set_yscale('log')
ax1.set_xlim(1e6*NP.min(intrfrmtr.lags)-1.0, 1e6*NP.max(intrfrmtr.lags)+1.0)
ax1.set_ylim(0.5*NP.min(NP.abs(intrfrmtr.skyvis_lag)),2.0*NP.max(NP.abs(intrfrmtr.skyvis_lag)))
l1, = ax1.plot([], [], 'g+', markersize=10)
ax1.tick_params(which='major', length=12, labelsize=18)
ax1.tick_params(which='minor', length=6)
# ax2 = fig.add_subplot(212)
# ax2.set_xlim(NP.min(skymod.catalog.location[:,0]), NP.max(skymod.catalog.location[:,0]))
# ax2.set_ylim(NP.min(skymod.catalog.location[:,1])-5.0, NP.max(skymod.catalog.location[:,1])+5.0)

ax2 = fig.add_subplot(212, projection='hammer')
ra_deg = skymod.catalog.location[:,0]
neg_ra = skymod.catalog.location[:,0] > 180.0
ra_deg[neg_ra] = ra_deg[neg_ra] - 360.0

ax2.set_xlabel(r'$\alpha$ [degrees]', fontsize=18)
ax2.set_ylabel(r'$\delta$ [degrees]', fontsize=18)
ax2.set_title('Sky Model', fontsize=18, weight='semibold')
# ax2.text(-2.0, -2.0, 'Sky Model', fontsize=18, va='bottom')
ax2.grid(True)
ax2.tick_params(which='major', length=12, labelsize=18)
ax2.tick_params(which='minor', length=6)

# l2init, = ax2.plot(skymod.catalog.location[:,0], skymod.catalog.location[:,1], 'k.', markersize=1)
l2init, = ax2.plot(NP.radians(ra_deg), NP.radians(skymod.catalog.location[:,1]), 'k.', markersize=1)

l2, = ax2.plot([], [], 'g+', markersize=3)

txt1 = ax1.text(0.05, 0.9, '', transform=ax1.transAxes, fontsize=18)
txt2 = ax2.text(0.25, 0.8, '', transform=ax2.transAxes, fontsize=18)

# def init():
#     l1.set_xdata([])
#     l1.set_ydata([])
#     l2.set_xdata(skymod.catalog.location[:,0])
#     l2.set_ydata(skymod.catalog.location[:,1])
#     l2.set_marker('.')
#     txt1.set_text('')
#     txt2.set_text('')
#     return l1, l2, txt1, txt2

def update(i, interferometer, line1, line2, t1, t2):
    line1.set_xdata(1e6 * interferometer.lags)
    line1.set_ydata(NP.abs(interferometer.vis_lag[i,:]))
    # line2.set_xdata(skymod.catalog.location[NP.asarray(interferometer.obs_catalog_indices[i]),0])
    # line2.set_ydata(skymod.catalog.location[NP.asarray(interferometer.obs_catalog_indices[i]),1])

    line2.set_xdata(NP.radians(ra_deg[NP.asarray(interferometer.obs_catalog_indices[i])]))
    line2.set_ydata(NP.radians(skymod.catalog.location[NP.asarray(interferometer.obs_catalog_indices[i]),1]))

    label_str = r' $\alpha$ = {0:+.3f} deg, $\delta$ = {1:+.2f} deg'.format(float(interferometer.timestamp[i])-interferometer.pointing_center[i,0], interferometer.pointing_center[i,1])
    t1.set_text(label_str)
    # t2.set_text(label_str)
    t2.set_text('')

    return line1, line2, t1, t2

anim = MOV.FuncAnimation(fig, update, fargs=(intrfrmtr, l1, l2, txt1, txt2), frames=intrfrmtr.vis_lag.shape[0], interval=400, blit=False)
PLT.show()
anim.save('/Users/t_nithyanandan/Downloads/delay_spectrum_animation.gif', fps=2.5, writer='imagemagick')
# anim.save('/Users/t_nithyanandan/Downloads/delay_spectrum_animation.mp4', fps=2.5, writer='ffmpeg')

# fig = PLT.figure(num=0, figsize=(14,12))
# ax1 = PLT.subplot2grid((2,1), (0,0))
# ax2 = PLT.subplot2grid((2,1), (1,0))

# ax1.set_title('Delay Spectrum', fontsize=12)
# ax2.set_title('Sky Model', fontsize=12)

# ax1.set_xlim(1e6*NP.min(intrfrmtr.lags)-1.0, 1e6*NP.max(intrfrmtr.lags)+1.0)
# ax1.set_ylim(0.5*NP.min(NP.abs(intrfrmtr.skyvis_lag)),2.0*NP.max(NP.abs(intrfrmtr.skyvis_lag)))
# ax1.set_yscale('log')

# ax2.set_xlim(NP.min(skymod.catalog.location[:,0]), NP.max(skymod.catalog.location[:,0]))
# ax2.set_ylim(NP.min(skymod.catalog.location[:,1])-5.0, NP.max(skymod.catalog.location[:,1])+5.0)

# ax2.grid(True)

# ax1.set_xlabel(r'$\eta$ [$\mu$s]', fontsize=12)
# ax1.set_ylabel('Amplitude [Jy Hz]', fontsize=12)
# ax2.set_xlabel(r'$\alpha$ [degrees]', fontsize=12)
# ax2.set_ylabel(r'$\delta$ [degrees]', fontsize=12)
# txt1 = ax1.text(0.05, 0.9, '', transform=ax1.transAxes, fontsize=12)
# txt2 = ax2.text(0.05, 0.9, '', transform=ax1.transAxes, fontsize=12)

# l1, = ax1.plot([], [], 'k+', markersize=10)
# l2init, = ax2.plot(skymod.catalog.location[:,0], skymod.catalog.location[:,1], 'k.', markersize=1)
# l2, = ax2.plot([], [], 'g+', markersize=2)

# # i = -1

# def init():
#     global skymod
#     l1.set_xdata([])
#     l1.set_ydata([])
#     l2.set_xdata(skymod.catalog.location[:,0])
#     l2.set_ydata(skymod.catalog.location[:,1])
#     l2.set_marker('.')
#     txt1.set_text('')
#     return l1, l2, txt1

# def updateData(i):
#     global intrfrmtr
#     # global i

#     # l1.set_xdata(1e6 * intrfrmtr.lags)
#     # l1.set_xdata(NP.abs(intrfrmtr.vis_lag[i,:]))
#     l1.set_data(1e6 * intrfrmtr.lags, NP.abs(intrfrmtr.vis_lag[i,:]))
#     l1.set_color('g')
#     # l2.set_xdata(skymod.catalog.location[NP.asarray(intrfrmtr.obs_catalog_indices[i]),0])
#     # l2.set_ydata(skymod.catalog.location[NP.asarray(intrfrmtr.obs_catalog_indices[i]),1])
#     l2.set_data(skymod.catalog.location[NP.asarray(intrfrmtr.obs_catalog_indices[i]),0], skymod.catalog.location[NP.asarray(intrfrmtr.obs_catalog_indices[i]),1])
#     l2.set_color('g')
#     label_str = 'RA = {0:+.3f} deg, DEC = {1:+.2f} deg'.format(float(intrfrmtr.timestamp[i])-intrfrmtr.pointing_center[i,0], intrfrmtr.pointing_center[i,1])    
#     txt1.set_text(label_str)
#     txt2.set_text(label_str)

#     return l1, l2, txt1

# anim = MOV.FuncAnimation(fig, updateData, blit=False, frames=intrfrmtr.vis_lag.shape[0], interval=400, repeat=True, init_func=init)
