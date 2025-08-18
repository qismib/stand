import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')


class BasicPlotter():

    tlist = []
    res = []
    ax_label = ()

    _num_plots = 1

    labels = None
    file_name = ""
    title = ""
    

    def __init__(self, times, results, ax_label):
        self.tlist = times
        self.res = results
        self._num_plots = len(self.res)
        self.ax_label = ax_label


    def save_plot(self):
        """
        save_plot
        """
        fig, ax = plt.subplots()
        ax.set_title(rf"{self.title}")
        
        if self.labels is not None and len(self.labels) == self._num_plots:
            for i in range (0, self._num_plots): ax.plot(self.tlist, self.res[i], label=self.labels[i])
            fig.legend(loc=1)
        else:
            for i in range (0, self._num_plots): ax.plot(self.tlist, self.res[i])

        ax.set_xlabel(self.ax_label[0])
        ax.set_ylabel(self.ax_label[1])

        plt.grid() # TODO: always on
        plt.tight_layout()
        plt.savefig(self.file_name)
        fig.clear()
        plt.close(fig)

