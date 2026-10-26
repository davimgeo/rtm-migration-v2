import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation
from scipy.ndimage import uniform_filter1d


def load(path: str, height: int, weight: int) -> np.ndarray:
  return np.fromfile(
    path, dtype=np.float32, count=height*weight
  ).reshape([height, weight], order='C')

def load_fortran(path: str, height: int, weight: int) -> np.ndarray:
  return np.fromfile(
    path, dtype=np.float32, count=height*weight
  ).reshape([height, weight], order='F')

def plot_model(model) -> None:
  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 8))

  img = ax.imshow(model, aspect="auto", cmap="jet")
  plt.colorbar(img, ax=ax, label="VP [m/s]")

  ax.set_title("Inverted Model", fontsize=13)
  ax.set_xlabel("Distance [m]", fontsize=13)
  ax.set_ylabel("Depth [m]", fontsize=13)

  plt.tight_layout()
  plt.show()

def plot_model_and_real(model, perc=99) -> None:
  real_model = load_fortran("data/FWI/marmousi_141z_681x_25dxdz.bin", 141, 681)

  fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))

  img = ax[0].imshow(model, aspect="auto", cmap="jet")
  plt.colorbar(img, ax=ax[0], label="VP [m/s]")

  img = ax[1].imshow(real_model, aspect="auto", cmap="jet")
  plt.colorbar(img, ax=ax[1], label="VP [m/s]")

  ax[0].set_title("Inverted Model", fontsize=13)
  ax[1].set_title("Real Model", fontsize=13)

  for ax in ax: 
    ax.set_xlabel("Distance [m]", fontsize=13)
    ax.set_ylabel("Depth [m]", fontsize=13)

  plt.tight_layout()
  plt.show()

def animate_iterations(
  iterations,
  shot=42,
  reference_trace=90,
  nt=12001,
  nrec=341,
  dt=0.001,
  nz=141,
  nx=681,
  tmin=1000,
  tmax=6500,
  receiver_spacing=50,
  interval=1000,
  save_path=None
):
  iterations = list(iterations)
  print(iterations)

  if not iterations:
    raise ValueError("iterations cannot be empty")

  dobs = load(
    DOBS_PATH +
    f"seismogram_12001x341_shot{shot}.bin",
    nt, nrec
  )

  traces = np.linspace(10, nrec - 11, 7, dtype=int)
  arg_min = np.argmin(np.abs(traces - reference_trace))
  traces[arg_min] = reference_trace

  time = np.arange(tmin, tmax) * dt
  window = max(1, int(0.5 / dt))

  models = [
    load(PATH + f"m_{iteration}.bin", nz, nx)
    for iteration in iterations
  ]

  vmin = min(np.percentile(model, 1) for model in models)
  vmax = max(np.percentile(model, 99) for model in models)

  def get_dcalc(iteration):
    return load(
      PATH + f"data/dcalc_{iteration}.bin",
      nt, nrec
    )

  def get_agc(dcalc):
    energy = 0.5 * (
      dcalc.astype(np.float64)**2
      + dobs.astype(np.float64)**2
    )

    rms = np.sqrt(
      uniform_filter1d(energy, size=window, axis=0)
    )

    floor = 0.02 * np.max(rms, axis=0, keepdims=True)
    scale = np.maximum(rms, floor + 1e-30)

    dcalc_agc = dcalc / scale
    dobs_agc = dobs / scale

    clip = np.percentile(
      np.abs(np.concatenate((
        dcalc_agc.ravel(),
        dobs_agc.ravel()
      ))),
      99
    )

    return dcalc_agc, dobs_agc, max(clip, 1e-10)

  dcalc = get_dcalc(iterations[0])
  dcalc_agc, dobs_agc, clip = get_agc(dcalc)

  fig = plt.figure(figsize=(20, 11), constrained_layout=True)
  gs = fig.add_gridspec(2, 7, height_ratios=[1, 1.2])

  axs = []

  for i in range(7):
    ax = fig.add_subplot(
      gs[0, i],
      sharey=axs[0] if axs else None
    )
    axs.append(ax)

  ax_model = fig.add_subplot(gs[1, :3])
  ax_calc = fig.add_subplot(gs[1, 3:5])
  ax_obs = fig.add_subplot(gs[1, 5:7])

  lines_calc = []
  lines_obs = []

  for i, irec in enumerate(traces):
    line_calc, = axs[i].plot(
      np.zeros_like(time), time,
      label="dcalc", linewidth=1
    )

    line_obs, = axs[i].plot(
      np.zeros_like(time), time,
      label="dobs", linewidth=1
    )

    lines_calc.append(line_calc)
    lines_obs.append(line_obs)

    offset = (irec - reference_trace) * receiver_spacing

    axs[i].set_title(
      f"Receiver {irec}\nOffset: {offset:+.0f} m",
      fontsize=10
    )

    axs[i].set_xlim(-1.1, 1.1)
    axs[i].set_ylim(time[-1], time[0])
    axs[i].grid(alpha=0.3)

    if i > 0:
      axs[i].tick_params(labelleft=False)

  axs[0].legend()
  axs[0].set_ylabel("Time (s)")
  axs[3].set_xlabel("Normalized Amplitude")

  img_model = ax_model.imshow(
    models[0],
    aspect="auto",
    cmap="jet",
    vmin=vmin,
    vmax=vmax
  )

  ax_model.set_title("Velocity Model")
  ax_model.set_xlabel("Distance [grid cells]")
  ax_model.set_ylabel("Depth [grid cells]")

  fig.colorbar(
    img_model, ax=ax_model,
    label="VP [m/s]"
  )

  extent = [
    0,
    (nrec - 1) * receiver_spacing,
    nt * dt,
    0
  ]

  img_calc = ax_calc.imshow(
    dcalc_agc,
    cmap="Greys",
    aspect="auto",
    vmin=-clip,
    vmax=clip,
    extent=extent
  )

  img_obs = ax_obs.imshow(
    dobs_agc,
    cmap="Greys",
    aspect="auto",
    vmin=-clip,
    vmax=clip,
    extent=extent
  )

  ax_calc.set_title("dcalc - AGC")
  ax_obs.set_title("dobs - AGC")

  for ax in [ax_calc, ax_obs]:
    ax.set_xlabel("Receiver Position (m)")

  ax_calc.set_ylabel("Time (s)")

  def update(frame):
    iteration = iterations[frame]

    dcalc = get_dcalc(iteration)

    for i, irec in enumerate(traces):
      trace_s = dcalc[tmin:tmax, irec]
      trace_o = dobs[tmin:tmax, irec]

      scale = max(
        np.max(np.abs(trace_s)),
        np.max(np.abs(trace_o))
      )

      if scale > 0:
        trace_s = trace_s / scale
        trace_o = trace_o / scale

      lines_calc[i].set_xdata(trace_s)
      lines_obs[i].set_xdata(trace_o)

    img_model.set_data(models[frame])

    dcalc_agc, dobs_agc, clip = get_agc(dcalc)

    img_calc.set_data(dcalc_agc)
    img_obs.set_data(dobs_agc)

    img_calc.set_clim(-clip, clip)
    img_obs.set_clim(-clip, clip)

    fig.suptitle(
      f"FWI - Iteration {iteration}",
      fontsize=17
    )

    return (
      lines_calc + lines_obs
      + [img_model, img_calc, img_obs]
    )

  update(0)

  ani = FuncAnimation(
    fig,
    update,
    frames=len(iterations),
    interval=interval,
    repeat=True,
    blit=False,
    cache_frame_data=False
  )

  if save_path is not None:
    writer = "pillow" if save_path.endswith(".gif") else "ffmpeg"

    ani.save(
      save_path,
      writer=writer,
      fps=1000 / interval
    )

  plt.show()

  return ani

def plot_agc(
  dobs,
  dcalc,
  dt=1e-3,
  receiver_spacing=50,
  window_time=0.5,
  clip_percentile=99
):
  nt, nrec = dobs.shape

  window = max(1, int(window_time / dt))

  energy = 0.5 * (
    dobs.astype(np.float64)**2
    + dcalc.astype(np.float64)**2
  )

  rms = np.sqrt(
    uniform_filter1d(energy, size=window, axis=0)
  )

  floor = 0.02 * np.max(rms, axis=0, keepdims=True)
  scale = np.maximum(rms, floor + 1e-30)

  dobs_agc = dobs / scale
  dcalc_agc = dcalc / scale

  clip = np.percentile(
    np.abs(np.concatenate((
      dobs_agc.ravel(),
      dcalc_agc.ravel()
    ))),
    clip_percentile
  )

  clip = max(clip, 1e-10)

  extent = [
    0,
    (nrec - 1) * receiver_spacing,
    (nt - 1) * dt,
    0
  ]

  fig, axs = plt.subplots(
    1, 2, figsize=(16, 8),
    sharey=True,
    constrained_layout=True
  )

  for ax, data, title in zip(
    axs,
    [dobs_agc, dcalc_agc],
    ["dobs - AGC", "dcalc - AGC"]
  ):
    ax.imshow(
      data,
      cmap="Greys",
      aspect="auto",
      vmin=-clip,
      vmax=clip,
      extent=extent
    )

    ax.set_title(title, fontsize=13)
    ax.set_xlabel("Receiver Position [m]")
    ax.set_ylabel("Time [s]")

  plt.show()

  return dobs_agc, dcalc_agc

PATH = "data/FWI/AWI_Test/"
DOBS_PATH = "data/FWI/dobs_12s_10hz_ng/"

MIN_ITERATION = 1
MAX_ITERATION = 21

#ani = animate_traces(
#  iteration=1,
#  shot=25,
#  tmin=1000,
#  tmax=12000,
#  t0=0.5,
#  receiver_step=1,
#  interval=100
#)

ani = animate_iterations(iterations=range(MIN_ITERATION, MAX_ITERATION + 1))

m_gradient = load("data/FWI/m0_gradient.bin", 141, 681)
m_current = load(PATH + f"m_{MAX_ITERATION}.bin", 141, 681)

#wiener = load("wiener.bin", 12001, 341)

#plot_model_and_real(m_current)
plot_model(m_current)

import numpy as np
import matplotlib.pyplot as plt

gradient = load("nabla_chi_141x681.bin", 141, 681)

# Your current velocity model in m/s
vp = m_gradient

nz, nx = gradient.shape

dz = 25.0
dx = 25.0
f_ref = 10.0

L = 0.3 * vp / f_ref

smooth_grad = np.zeros_like(gradient)

for i in range(nz):
  for j in range(nx):
    length = L[i, j]

    sigma_z = length / dz
    sigma_x = length / dx

    rz = int(np.ceil(3 * sigma_z))
    rx = int(np.ceil(3 * sigma_x))

    i0 = max(0, i - rz)
    i1 = min(nz, i + rz + 1)

    j0 = max(0, j - rx)
    j1 = min(nx, j + rx + 1)

    z = (np.arange(i0, i1) - i) * dz
    x = (np.arange(j0, j1) - j) * dx

    G = np.exp(
      -0.5 * (
        (z[:, None] / length)**2 +
        (x[None, :] / length)**2
      )
    )

    window = gradient[i0:i1, j0:j1]

    smooth_grad[i, j] = np.sum(G * window) / np.sum(G)

fig, axs = plt.subplots(1, 3, figsize=(14, 5))

axs[0].imshow(gradient, cmap="Greys", aspect="auto")
axs[0].set_title("Original Gradient")

axs[1].imshow(smooth_grad, cmap="Greys", aspect="auto")
axs[1].set_title("Gaussian-Smoothed Gradient")

axs[2].imshow(gradient - smooth_grad, cmap="Greys", aspect="auto")
axs[2].set_title("Difference")

plt.tight_layout()
plt.show()
