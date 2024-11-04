class Foo():
    def __init__(self, a, b):
        self.parameters = {key: value for key, value in locals().items() if not key.startswith('__') and key != 'self'}
        self.a = a
        self.b = b

    def bar(self):
        print("bar", self.a)

    def zop(self):
        print("zip", self.b)

foo = Foo(1,2)
print(foo.parameters)

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

class SliderPlot:
    def __init__(self, title='Slider Example'):
        self.fig, self.ax = plt.subplots()
        self.fig.subplots_adjust(left=0.2, bottom=0.25)  # Adjust to make room for sliders
        self.sliders = []
        self.slider_y_position = 0.1  # Initial y position for the first slider
        self.title = title
        plt.title(self.title)

    def add_slider(self, label, slider_range, initial_value):
        # Create a new slider
        ax_slider = self.fig.add_axes([0.1, self.slider_y_position, 0.65, 0.03])  # Position and size
        slider = Slider(ax_slider, label, slider_range[0], slider_range[1], valinit=initial_value)

        # Update the y position for the next slider
        self.slider_y_position += 0.05  # Move down for the next slider
        self.sliders.append(slider)

        return slider  # Return the slider object if needed

    def show(self):
        plt.show()

# Example usage
if __name__ == "__main__":
    slider_plot = SliderPlot()
    slider_plot.add_slider('Slider 1', (0, 10), 5)
    slider_plot.add_slider('Slider 2', (0, 100), 50)
    slider_plot.add_slider('Slider 3', (0, 1), 0.5)
    slider_plot.show()