import numpy as np 
import matplotlib.pyplot as plt

def plot_seismogram(
    seismogram: np.ndarray, 
    dt: float, 
    offset: int, 
    perc=99
) -> None:

  nt, nrec = seismogram.shape
  
  tloc = np.linspace(0, nt - 1, 11, dtype=int)
  tlab = np.around(tloc * dt, decimals=1)

  xloc = np.linspace(0, nrec - 1, 9)
  xlab = np.array(offset * 10 * xloc, dtype=int)

  scale_min = np.percentile(seismogram, 100 - perc)
  scale_max = np.percentile(seismogram, perc)

  fig, ax = plt.subplots(figsize=(10, 8))

  img = ax.imshow(seismogram, aspect="auto", cmap="Greys",
                    vmin=scale_min, vmax=scale_max)

  ax.set_yticks(tloc)
  ax.set_yticklabels(tlab)

  ax.set_xticks(xloc)
  ax.set_xticklabels(xlab)

  ax.set_xlabel("Offset (m)", fontsize=13)
  ax.set_ylabel("TWT (s)", fontsize=13)

  plt.show()

def plot_seismogram_elastic(
    calc_p: np.ndarray, 
    vx: np.ndarray,
    vz: np.ndarray,
    dt: float, 
    offset: int, 
    perc=99
) -> None:

  nt, nrec = calc_p.shape
  
  tloc = np.linspace(0, nt - 1, 7, dtype=int)
  tlab = np.around(tloc * dt, decimals=1)

  xloc = np.linspace(0, nrec - 1, 7)
  xlab = np.array(offset * 10 * xloc, dtype=int)

  fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(15, 8))

  img = ax[0].imshow(calc_p, aspect="auto", cmap="Greys")
  ax[0].set_title("Pressure", fontsize=15)
  plt.colorbar(img, ax=ax[0])

  img1 = ax[1].imshow(vx, aspect="auto", cmap="Greys")
  ax[1].set_title("Vx", fontsize=15)
  plt.colorbar(img1, ax=ax[1])

  img2 = ax[2].imshow(vz, aspect="auto", cmap="Greys")
  ax[2].set_title("Vz", fontsize=15)
  plt.colorbar(img2, ax=ax[2])

  for axs in ax:
    axs.set_yticks(tloc)
    axs.set_yticklabels(tlab)

    axs.set_xticks(xloc)
    axs.set_xticklabels(xlab)

    axs.set_xlabel("Offset (m)", fontsize=13)
    axs.set_ylabel("TWT (s)", fontsize=13)
  
  plt.tight_layout()
  plt.show()

def plot_model(model, perc=99) -> None:
  _, ax = plt.subplots(figsize=(12, 8))

  vmin = np.percentile(model, 100 - perc)
  vmax = np.percentile(model, perc)

  img = ax.imshow(
    model, aspect="auto", cmap="jet", 
    vmin=vmin, vmax=vmax
  )

  plt.colorbar(img, ax=ax, label="VP [m/s]")

  ax.set_xlabel("Distance [m]", fontsize=13)
  ax.set_ylabel("Depth [m]", fontsize=13)
  ax.set_title("Velocity Model", fontsize=13)

  plt.tight_layout()
  plt.show()

def plot_model_geometry(
  model: np.ndarray, 
  nb: int,
  dh: int, 
  recx: np.ndarray, 
  recz: np.ndarray, 
  srcx: np.ndarray, 
  srcz: np.ndarray,
  perc=99
) -> None:
  nzz, nxx = model.shape

  nz = nzz - 2*nb
  nx = nxx - 2*nb

  model = model[nb:nb + nz, nb:nb + nx]

  _, ax = plt.subplots(figsize=(12, 5))

  vmin = np.percentile(model, 100 - perc)
  vmax = np.percentile(model, perc)

  xloc = np.linspace(0, nx - 1, 11, dtype=int)
  xlab = np.array(xloc * dh, dtype=int)

  zloc = np.linspace(0, nz - 1, 7, dtype=int)
  zlab = np.array(zloc * dh, dtype=int)

  img = ax.imshow(
    model, aspect="auto", cmap="jet", 
    vmin=vmin, vmax=vmax
  )

  ax.plot(recx, recz, 'gv', label="Receivers", markersize=12)
  ax.plot(srcx, srcz, 'r*', markersize=12, label="Source")

  plt.colorbar(img, ax=ax, label="VP [m/s]")

  ax.set_xlabel("Distance [m]", fontsize=13)
  ax.set_ylabel("Depth [m]", fontsize=13)
  ax.set_title("Velocity Model", fontsize=13)

  plt.tight_layout()
  plt.show()

def plot1d(arr: np.ndarray) -> None:
  _, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 5))

  ax.plot(arr)
  
  plt.tight_layout()
  plt.show()

def plot2d(arr: np.ndarray) -> None:
  _, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 5))

  ax.imshow(arr, aspect="auto", cmap="jet")
  
  plt.tight_layout()
  plt.show()
