class Spacecraft:
    def __init__(self, fuel_tanks, oxidizer_tanks):
        self.fuel_tanks = fuel_tanks
        self.oxidizer_tanks = oxidizer_tanks

    def total_mass(self):
        total = sum(tank.weight for tank in self.fuel_tanks + self.oxidizer_tanks)
        return total