import yaml

class PlotConfig:
    """Configuration for plotting figures.

    Attributes:
        figsize: size of figure
        facecolor: colour of background
        grid: toggles grid visibility
        grid_color: color of grid lines
        xtick_width: frequency of xticks
        xtick_size: size of ticks on x axis
        ytick_width: frequency of yticks
        ytick_size: size of ticks on y axis
        linewidth: pixel width of line on plot
        fontsize: size of text labels
        xtick_btm: toggles tick positions
        ytick_left: toggle ytick position
        show: toggles showing plot on screen
        save: toggles saving figure to file

    """

    def __init__(self, input_file):
        with open(input_file, "r") as ymlfile:
            inputs = yaml.load(ymlfile, Loader=yaml.FullLoader)

        self.figsize = inputs["PLOT"]["FIG_SIZE"]
        self.facecolor = inputs["PLOT"]["FACECOLOR"]
        self.grid = inputs["PLOT"]["GRID"]
        self.grid_color = inputs["PLOT"]["GRIDCOLOR"]
        self.xtick_width = inputs["PLOT"]["XTICK_WIDTH"]
        self.xtick_size = inputs["PLOT"]["XTICK_SIZE"]
        self.ytick_width = inputs["PLOT"]["YTICK_WIDTH"]
        self.ytick_size = inputs["PLOT"]["YTICK_SIZE"]
        self.linewidth = inputs["PLOT"]["LINEWIDTH"]
        self.fontsize = inputs["PLOT"]["FONTSIZE"]
        self.xtick_btm = inputs["PLOT"]["XTICK_BTM"]
        self.ytick_left = inputs["PLOT"]["YTICK_LEFT"]
        self.save = inputs["PLOT"]["SAVE"] 