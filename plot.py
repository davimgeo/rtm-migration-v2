import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def load(path: str, height: int, weight: int) -> np.ndarray:
  return np.fromfile(
    path, dtype=np.float32, count=height*weight
  ).reshape([height, weight], order='C')

nt = 8001
nrec = 341

dcalc = load("data/dcalc.bin", nt, nrec)
dobs = load("data/dobs.bin", nt, nrec)
diff = dobs - dcalc

plt.plot(diff[:, 75])
plt.plot(diff[:, 300])
plt.show()

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 8))

ax[0].imshow(dcalc, aspect="auto", cmap="Greys")
trace_dcalc, = ax[0].plot(
  np.zeros(nt),
  np.arange(nt),
  '--'
)

ax[1].imshow(dobs, aspect="auto", cmap="Greys")
trace_dobs, = ax[1].plot(
  np.zeros(nt),
  np.arange(nt),
  '--'
)

line_dcalc, = ax[2].plot([], [], label="dcalc")
line_dobs, = ax[2].plot([], [], label="dobs")
line_diff, = ax[2].plot([], [], label="diff")

ax[2].set_xlim(
  min(dcalc.min(), dobs.min(), diff.min()),
  max(dcalc.max(), dobs.max(), diff.max())
)
ax[2].set_ylim(nt - 1, 0)

ax[2].legend()

title = fig.suptitle("Trace 0")

def update(trace):
  trace_dcalc.set_xdata(np.ones(nt) * trace)
  trace_dobs.set_xdata(np.ones(nt) * trace)

  line_dcalc.set_data(dcalc[:, trace], np.arange(nt))
  line_dobs.set_data(dobs[:, trace], np.arange(nt))
  line_diff.set_data(diff[:, trace], np.arange(nt))

  title.set_text(f"Trace {trace}")

  return (
    trace_dcalc,
    trace_dobs,
    line_dcalc,
    line_dobs,
    line_diff
  )

ani = FuncAnimation(
  fig,
  update,
  frames=range(nrec),
  interval=50,
  repeat=True
)

plt.tight_layout()
plt.show()
