import sys
import pandas
import matplotlib.pyplot as plt

logfile = sys.argv[1]

pddata = pandas.read_csv(logfile, delim_whitespace=True, comment="#")


fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(20,5))

ax1.semilogy(pddata["iter"], pddata["res_c"], label="Residual")
ax1.set_xlabel("Iteration")
ax1.set_ylabel("Residual")
ax1.legend()

ax2.plot(pddata["iter"], pddata["min_c"], label="Min C")
ax2.plot(pddata["iter"], pddata["max_c"], label="Max C")
ax2.set_xlabel("Iteration")
ax2.set_ylabel("c")
ax2.legend()

ax3.plot(pddata["iter"], pddata["grad_c_mag"], label="||grad(c)||")
ax3.set_xlabel("Iteration")
ax3.set_ylabel("||grad(c)||")
ax3.legend()

plt.show()
