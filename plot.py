import numpy as np 
import matplotlib.pyplot as plt 

def load(path: str, height: int, weight: int) -> np.ndarray:
  return np.fromfile(
    path, dtype=np.float32, count=height*weight
  ).reshape([height, weight], order='C')

dcalc = load("data/dcalc.bin", 4001, 114)
dobs = load("data/dobs.bin", 4001, 114)
diff = dobs - dcalc

trace = 15

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 8))

ax[0].imshow(dcalc, aspect="auto", cmap="Greys")
ax[0].plot(np.ones(4001)*trace, np.arange(4001), '--')

ax[1].imshow(dobs, aspect="auto", cmap="Greys")
ax[1].plot(np.ones(4001)*trace, np.arange(4001), '--')

ax[2].plot(dcalc[:, trace], label="dcalc")
ax[2].plot(dobs[:, trace], label="dobs")
plt.plot(diff[:, trace], label="diff")

plt.legend()
plt.tight_layout()
plt.show()

