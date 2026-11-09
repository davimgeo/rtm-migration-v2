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

def animate_iterations(
  iterations,
  shot=25,
  reference_trace=90,
  nt=6001,
  nrec=101,
  dt=0.001,
  nz=101,
  nx=101,
  tmin=500,
  tmax=3000,
  receiver_spacing=10,
  interval=1000,
  save_path=None
):
  iterations = list(iterations)
  print(iterations)

  if not iterations:
    raise ValueError("iterations cannot be empty")

  dobs = load(
    DOBS_PATH +
    f"seismogram_6001x101_shot{shot}.bin",
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

PATH = "data/FWI/AWI_Circle_1700_2/"
DOBS_PATH = "data/dobs/"

MIN_ITERATION = 1
MAX_ITERATION = 5

ani = animate_iterations(iterations=range(MIN_ITERATION, MAX_ITERATION + 1))

model = load(PATH + f"m_{MAX_ITERATION}.bin", 101, 101)
plot_model(model)

phi_history = np.loadtxt(PATH + "phi_history.txt", dtype=np.float32)
max_phi = np.max(phi_history)
phi_history /= max_phi

plt.plot(phi_history, '-o', label="gaussian")

plt.grid(True)
plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt

wiener = load(PATH + "data/wiener_1.bin", 6001, 101)

rec = 40
dt = 1e-3
sigma = 0.30

w1 = np.fft.fftshift(wiener[:, rec].copy())
w1 /= max(np.max(np.abs(w1)), 1e-12)

nt = len(w1)
lags = np.arange(nt) - nt // 2
time = lags * dt

P = np.exp(-(time**2) / sigma**2)

w_penalized = P * w1

fig, ax = plt.subplots(
    nrows=1, ncols=3,
    figsize=(18, 5),
    sharex=True
)

ax[0].plot(time, w1)
ax[0].set_title("Original Wiener (Normalized)")

ax[1].plot(time, P)
ax[1].set_title("Penalty: 1 - Gaussian")

ax[2].plot(time, w1, "--", alpha=0.4, label="Original")
ax[2].plot(time, w_penalized, label="Penalized")
ax[2].set_title("Wiener × Penalty")
ax[2].legend()

for a in ax:
    a.set_xlabel("Lag [s]")
    a.set_ylabel("Amplitude")
    a.grid(alpha=0.3)
    a.axvline(0, color="red", linestyle="--", alpha=0.5)
    a.set_xlim(-2, 2)

plt.tight_layout()
plt.show()

model1 = load(PATH + f"m_{MAX_ITERATION}.bin", 141, 681)
model2 = load("data/FWI/L2_Good_Test2/" + f"m_{MAX_ITERATION}.bin", 141, 681)

diff = model1 - model2

fig, ax = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

vmin = min(model1.min(), model2.min())
vmax = max(model1.max(), model2.max())

dmax = np.max(np.abs(diff))

im0 = ax[0].imshow(model1, cmap="jet", aspect="auto")
ax[0].set_title("L2 Good 1")
fig.colorbar(im0, ax=ax[0])

im1 = ax[1].imshow(model2, cmap="jet", aspect="auto")
ax[1].set_title("L2 Good 2")
fig.colorbar(im1, ax=ax[1])

im2 = ax[2].imshow(diff, cmap="seismic", vmin=-dmax, vmax=dmax, aspect="auto")
ax[2].set_title("Difference (Model 1 - Model 2)")
fig.colorbar(im2, ax=ax[2])

plt.show()

