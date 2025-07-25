import numpy as np
class pulse:
    def __init__(self, pulse_width, pulse_peak_time):
        self.pulse_width = pulse_width
        self.pulse_peak_time = pulse_peak_time

    def get_pulse_width(self):
        return self.pulse_width

    def get_pulse_peak_time(self):
        return self.pulse_peak_time

    def sample_launch_time(self):
        return np.random.normal(loc=self.pulse_peak_time, scale=self.pulse_width / 2.355)


