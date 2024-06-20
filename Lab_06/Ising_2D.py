"""
2D Ising Model: I set up a grid of NxN spins, White cell is Up and Black is Down
In Ising Model ther's only the interaction between nearest neighbor spins.
I compute it using Metropolis algorithm to sample Boltzmann weight 
"""
from manim import*
from colour import Color
from random import uniform

config.frame_rate = 60

class Background(VGroup): # Spin-Grid inherit from VGroup
    def __init__(self, bkg_len = config.frame_height-0.5, n_quad = 40, J = 1, T = 2.269 , color = GRAY):
        super().__init__()
        self.n_quad = n_quad # number of spins per row/column
        self.side_len = bkg_len / n_quad # hight = width of a square
        self.J = J # Spins interaction constant
        self.Beta = 1/T
        self.color = color # color of the grid
        self.set_grid() # set the grid

    def set_grid(self): # method that creates the grid of squares
        for _ in range(self.n_quad):
            self.add(VGroup(*[Square( # adding every row to the self VGroup
                side_length = self.side_len,
                color = BLACK, # I set all the starting squares of the same color (same down spin, e.g. zero temperature)
                fill_opacity = 1,
                stroke_color = self.color,
                stroke_width = self.side_len)
                for _ in range(self.n_quad)]).arrange(RIGHT, buff=0)    
            ).arrange(DOWN, buff=0)

    def spin_up(self, row, col): # method that turns a spin to up value
        self[row][col].set(fill_color = WHITE) # I use set method because I want to change the inner color only (not the grid wedge)
    
    def spin_down(self, row, col): # method that turns a spin to down value
        self[row][col].set(fill_color = BLACK)
    
    def change_spin(self, row, col): # method that flip the spin in the given position
        if self[row][col].get_color() == Color(WHITE): # spin UP
            self.spin_down(row, col)
        elif self[row][col].get_color() == Color(BLACK): # spin DOWN
            self.spin_up(row, col)
        else:
            print("Error: color grid must be White or Black")
            exit(0)

    def SpinFlipper(self, N_updates=1): # method that turns a certain spin if Metropolis criteria is sutisfied
        for _ in range(N_updates): # we can choose the number of interactions per function call   
            row = int(uniform(0, self.n_quad - 1))  # Random row used to select the lattice cell to analyze
            col = int(uniform(0, self.n_quad - 1))  # Random column
            # computing the difference in energy flipping that spin
            if self.DeltaEnergy(row, col) <=0: # System energy reduced -> we accept the new configuration
                self.change_spin(row, col)
            else:
                r = uniform(0, 1)
                A = np.exp(-self.Beta*self.DeltaEnergy(row, col)) # Acceptance value, Boltzmann weight
                if r < A: # We accept the new configuration
                    self.change_spin(row, col)
        # analogue, shorter but it has to calculate more values every call:
            # if (min(1, np.exp(-self.Beta*self.DeltaEnergy(row, col))) <= uniform(0,1)): self.change_spin(row, col)
        return self

    def DeltaEnergy(self, row, col): # method that computes E_mu - E_nu, the difference in energy between the newer state and the current one
        # Computing the spin value in the given lattice coordinates
        if self[row][col].get_color() == Color(WHITE): # spin UP
            self.actual_spin = 1
        elif self[row][col].get_color() == Color(BLACK): # spin DOWN
            self.actual_spin = -1
        else:
            print("Error: color grid must be White or Black")
            exit(0)
        # Energy difference between a given state and another one with flipped row,col spin
        energy_difference = 0

        # Calc of the energy of the near neighbors spin only
        for i in range(row - 1, row + 2):
            for j in range(col - 1, col + 2):
                # Periodic boundary conditions
                neighbor_row = (i + self.n_quad) % self.n_quad
                neighbor_col = (j + self.n_quad) % self.n_quad

                # Computing energy difference
                if (neighbor_row, neighbor_col) != (row, col):
                    # Invert the spin value in the given lattice coordinates
                    if self[neighbor_row][neighbor_col].get_color() == Color(WHITE): # spin UP
                        self.neighbor_spin = 1
                    elif self[neighbor_row][neighbor_col].get_color() == Color(BLACK): # spin DOWN
                        self.neighbor_spin = -1
                    else:
                        print("Error: color grid must be White or Black")
                        exit(0)
                    energy_difference += 2 * self.J * self.neighbor_spin * self.actual_spin
        return energy_difference

    def initialization(self, start_coords = []): # method that set the spin at the initial given value, I have to put spin up coordinates like this: [(row_1, col_1), (row_2, col2), ...]
        for r, c in start_coords:
            self.spin_up(r, c)

    def random_initialization(self): # method that set a random configuration for the grid
        for row in range(self.n_quad):
            for col in range(self.n_quad):
                rand = uniform(0,1)
                if rand < 0.5:
                    self.change_spin(row, col)
        
class CallCountDTDecorator: # Decorator defined to set the correct value of dt, used in Manim "dt updaters", it counts how many time a function is called
    def __init__(self, inline_func):
        self.call_count = 0     # Counter
        self.inline_func = inline_func

    def __call__(self, dt, dt_calculate): # This set the initial dt (dt=0) to the real dt (1/fps)
        self.call_count = 0 if dt > 0 else self.call_count + 1
        if dt ==  0 and self.call_count > 1:
            self.call_count = 0
            return dt_calculate # it returns the real value of dt
        return dt # Althought it returns the dt value given by the dt updater  

@CallCountDTDecorator
def fix_dt(dt, dt_calculate):
    pass

class IsingModel(Scene):
    PAUSE = 0 # Updating time of the grid (must be > 1/fps)
    DURATION = 10 # Total time of the animation
    STEP_PAUSE = PAUSE # Integer multiple of "PAUSE"

    def construct(self):
        DT = 1 / self.camera.frame_rate # REAL dt
        spin_board = Background(T=0.02, n_quad = 100) # Saving the spins grid
        spin_board.random_initialization() # Initializing random spin up
        # spin_board.initialization(start_coords=[(20,20)])

        spin_board.self_time = 0 # Creating ad attribute used to trace the passing time
        
        def Spin_updater(mob, dt): # updater
            mob.self_time += fix_dt(dt, DT) # Accumulating proper time
            if mob.self_time > self.STEP_PAUSE: # I update the system every moment the proper time is grater than the pause
                mob.SpinFlipper(N_updates = 20) # Update the system
                self.STEP_PAUSE += self.PAUSE # Update the pause wall
             
        spin_board.add_updater(Spin_updater)
        self.add(spin_board)

        self.wait(self.DURATION)

        spin_board.suspend_updating()
        self.wait()