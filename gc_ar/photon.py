class Photon:
    def __init__(self, position, direction, stokes,energy):
        self.position = position
        self.direction = direction
        self.stokes = stokes
        self.energy = energy

    def move(self, step_size):
        self.position += step_size * self.direction

    def rotate_stokes(self, matrix):
        self.stokes = matrix @ self.stokes

    def rotate_direction(self, matrix):
        self.direction = matrix @ self.direction

    def decay(self,energy):
        self.energy = energy

