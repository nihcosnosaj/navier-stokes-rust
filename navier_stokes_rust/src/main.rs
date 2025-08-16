
use piston_window::*;

pub struct FluidGrid {
    // Number of cells in x and y direction.
    nx: usize,
    ny: usize,

    // Cell width/height
    dx: f64,

    // Pressure (stored at cell centers)
    // Size: nx * ny
    p: Vec<f64>,

    // Horizontal velocity (stored on vertical faces)
    // Size: (nx + 1) * ny
    u: Vec<f64>,

    // Vertical velocity (stored on horizontal faces)
    // Size: nx * (ny + 1)
    v: Vec<f64>,
}

use piston_window::{Context, G2d, line, Transformed};
impl FluidGrid {
    // A function to create a new, empty grid.
    pub fn new(nx: usize, ny: usize, dx: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            p: vec![0.0; nx * ny],
            u: vec![0.0; (nx + 1) * ny],
            v: vec![0.0; nx * (ny + 1)],
        }
    }

    // Convert 2D u-velocity index into 1D vector index.
    fn u_idx(&self, i: usize, j: usize) -> usize {
        j * (self.nx + 1) + i
    }

    // Convert 2D v-velocity index to 1D vector index
    fn v_idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    fn get_velocity(&self, x: f64, y: f64) -> (f64, f64) {
        // Clamp position to be within the grid
    let x = x.max(0.0).min((self.nx as f64) * self.dx);
    let y = y.max(0.0).min((self.ny as f64) * self.dx);

        // --- Interpolate U-velocity ---
        // Find bottom-left u-velocity index and interpolation weights
        let u_i = (x / self.dx).floor() as usize;
        let u_j = ((y / self.dx) - 0.5).floor() as usize;
        let u_tx = (x / self.dx) - u_i as f64;
        let u_ty = (y / self.dx - 0.5) - u_j as f64;

        // Get the four surrounding u-velocities
        let u00 = self.u[self.u_idx(u_i.min(self.nx-1), u_j.min(self.ny-1))];
        let u10 = self.u[self.u_idx((u_i + 1).min(self.nx-1), u_j.min(self.ny-1))];
        let u01 = self.u[self.u_idx(u_i.min(self.nx-1), (u_j + 1).min(self.ny-1))];
        let u11 = self.u[self.u_idx((u_i + 1).min(self.nx-1), (u_j + 1).min(self.ny-1))];

        // Perform bilinear interpolation for u
        let u_interp = u00 * (1.0 - u_tx) * (1.0 - u_ty)
                     + u10 * u_tx * (1.0 - u_ty)
                     + u01 * (1.0 - u_tx) * u_ty
                     + u11 * u_tx * u_ty;

        // --- Interpolate V-velocity ---
        // Find bottom-left v-velocity index and interpolation weights
        let v_i = ((x / self.dx) - 0.5).floor() as usize;
        let v_j = (y / self.dx).floor() as usize;
        let v_tx = (x / self.dx - 0.5) - v_i as f64;
        let v_ty = (y / self.dx) - v_j as f64;

        // Get the four surrounding v-velocities
        let v00 = self.v[self.v_idx(v_i.min(self.nx-1), v_j.min(self.ny-1))];
        let v10 = self.v[self.v_idx((v_i + 1).min(self.nx-1), v_j.min(self.ny-1))];
        let v01 = self.v[self.v_idx(v_i.min(self.nx-1), (v_j + 1).min(self.ny-1))];
        let v11 = self.v[self.v_idx((v_i + 1).min(self.nx-1), (v_j + 1).min(self.ny-1))];

        // Perform bilinear interpolation for v
        let v_interp = v00 * (1.0 - v_tx) * (1.0 - v_ty)
                     + v10 * v_tx * (1.0 - v_ty)
                     + v01 * (1.0 - v_tx) * v_ty
                     + v11 * v_tx * v_ty;

        (u_interp, v_interp)
    }


    // Step 1: Advection
    // Moves the velocity field along itself.
    fn advect(&mut self, dt: f64) {
        // Implementation of semi-lagrangian advection will go here.
        // We'll trace back and interpolate for each u and v componnet.
    let mut u_new = self.u.clone();
    let mut v_new = self.v.clone();

        // Advect u-velocity components (excluding boundaries for simplicity)
        for j in 0..self.ny {
            for i in 1..self.nx {
                // Find the real-world position (x,y) of this u-velocity component
                let x = i as f64 * self.dx;
                let y = (j as f64 + 0.5) * self.dx;

                // Get the velocity at this exact point
                let (u, v) = self.get_velocity(x, y);

                // Trace back in time to find the source position
                let x_prev = x - dt * u;
                let y_prev = y - dt * v;

                // Get the advected velocity by interpolating at the source position
                let (u_advected, _) = self.get_velocity(x_prev, y_prev);
                u_new[self.u_idx(i, j)] = u_advected;
            }
        }

        // Advect v-velocity components (excluding boundaries for simplicity)
        for j in 1..self.ny {
            for i in 0..self.nx {
                // Find the real-world position (x,y) of this v-velocity component
                let x = (i as f64 + 0.5) * self.dx;
                let y = j as f64 * self.dx;

                // Get the velocity at this exact point
                let (u, v) = self.get_velocity(x, y);

                // Trace back in time to find the source position
                let x_prev = x - dt * u;
                let y_prev = y - dt * v;

                // Get the advected velocity by interpolating at the source position
                let (_, v_advected) = self.get_velocity(x_prev, y_prev);
                v_new[self.v_idx(i, j)] = v_advected;
            }
        }

        // Update the grid's velocity fields
        self.u = u_new;
        self.v = v_new;
    }

    fn p_idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    // Step 2: Pressure Solve
    // Solves the Poisson equation to enforce incompressibility.
    fn solve_pressure(&mut self, dt: f64) {
        // Implementation of the Jacobi iteration for the pressure solve.
        // This is where we calculate divergence and iterate to find p.
        let dx = self.dx;
        // We can assume density is 1 for simplicity, as it scales the pressure
        let rho = 1.0; 

        // Part 1: Calculate the divergence of the velocity field.
        // This is the right-hand side (RHS) of our Poisson equation.
        let mut divergence = vec![0.0; self.nx * self.ny];
        for j in 0..self.ny {
            for i in 0..self.nx {
                let u_right = self.u[self.u_idx(i + 1, j)];
                let u_left  = self.u[self.u_idx(i, j)];
                let v_top   = self.v[self.v_idx(i, j + 1)];
                let v_bot   = self.v[self.v_idx(i, j)];

                let d = (u_right - u_left + v_top - v_bot) / dx;
                divergence[self.p_idx(i, j)] = d;
            }
        }

        // Part 2: Iteratively solve for pressure using the Jacobi method.
        // We repeat this loop to let the pressure values settle.
        let mut p_new = self.p.clone();
        let num_iterations = 50; // More iterations = more accuracy
        for _ in 0..num_iterations {
            for j in 1..self.ny - 1 { // We only solve for interior pressure points
                for i in 1..self.nx - 1 {
                    let p_right = self.p[self.p_idx(i + 1, j)];
                    let p_left  = self.p[self.p_idx(i - 1, j)];
                    let p_top   = self.p[self.p_idx(i, j + 1)];
                    let p_bot   = self.p[self.p_idx(i, j - 1)];

                    let d = divergence[self.p_idx(i, j)];

                    // This is the discretized Poisson equation rearranged for p_i,j
                    let p_updated = (p_right + p_left + p_top + p_bot - d * dx * dx) / 4.0;
                    p_new[self.p_idx(i, j)] = p_updated;
                }
            }
            // Update the pressure field for the next iteration
            self.p = p_new.clone();
        }     

    }

    // Step 3: Projection
    // Corrects the velocity field using the pressure gradient.
    fn project(&mut self, dt: f64) {
        // Implementation of final velocity update.
        // We subtract the pressure gradient from the advected velocity.
            let dx = self.dx;
        let rho = 1.0; // Must match the density used in solve_pressure
        let scale = dt / (rho * dx);

        // Update the u-velocity components
        // We iterate over the interior faces where u is defined
        for j in 1..self.ny - 1 {
            for i in 1..self.nx {
                // Get the pressure values in the cells to the left and right of the u-component
                let p_left = self.p[self.p_idx(i - 1, j)];
                let p_right = self.p[self.p_idx(i, j)];

                // Apply the correction
                let u_idx = self.u_idx(i, j);
                self.u[u_idx] -= scale * (p_right - p_left);
            }
        }

        // Update the v-velocity components
        // We iterate over the interior faces where v is defined
        for j in 1..self.ny {
            for i in 1..self.nx - 1 {
                // Get the pressure values in the cells below and above the v-component
                let p_bot = self.p[self.p_idx(i, j - 1)];
                let p_top = self.p[self.p_idx(i, j)];

                // Apply the correction
                let v_idx = self.v_idx(i, j);
                self.v[v_idx] -= scale * (p_top - p_bot);
            }
        }
   
    }

    fn set_boundaries(&mut self) {
        // --- Vertical Walls (Left and Right) ---
        // Set u-velocity to 0 on the left and right walls.
        let nx = self.nx;
        let ny = self.ny;
        for j in 0..ny {
            let left = self.u_idx(0, j);
            let right = self.u_idx(nx, j);
            self.u[left] = 0.0;
            self.u[right] = 0.0;
        }

        // --- Horizontal Walls (Top and Bottom) ---
        // Set v-velocity to 0 on the top and bottom walls.
        for i in 0..nx {
            let bottom = self.v_idx(i, 0);
            let top = self.v_idx(i, ny);
            self.v[bottom] = 0.0;
            self.v[top] = 0.0;
        }
    }

    pub fn run_step(&mut self, dt: f64) {
        self.advect(dt);
        self.solve_pressure(dt);
        self.project(dt);
        self.set_boundaries();
    }

    pub fn draw_velocities(&self, c: &Context, g: &mut G2d) {
        // Bright yellow, 3px thick
        let line = line::Line::new([1.0, 1.0, 0.0, 1.0], 3.0);
        let scale = 2.0; // Larger scaling factor for arrow length

        // Loop over a slightly smaller grid to avoid drawing on the boundaries
        for j in 1..self.ny - 1 {
            for i in 1..self.nx - 1 {
                // Calculate the center of the cell
                let x = (i as f64 + 0.5) * self.dx;
                let y = (j as f64 + 0.5) * self.dx;

                // Get the velocity at that point
                let (u, v) = self.get_velocity(x, y);

                // Calculate the start and end points of the arrow
                let x2 = x + u * scale;
                let y2 = y + v * scale;

                // Draw the line
                line.draw([x, y, x2, y2], &c.draw_state, c.transform, g);
            }
        }
    }
}



fn main() {
    let mut grid = FluidGrid::new(40, 40, 15.0);
    // TODO: set some initial fluid motion.
    let center_i = grid.nx / 2;
    let center_j = grid.ny / 2;
    let idx = grid.v_idx(center_i, center_j);
    grid.v[idx] = 100.0;

    let mut window: PistonWindow = WindowSettings::new("Fluid Sim", [600, 600]).exit_on_esc(true).build().unwrap();

    while let Some(event) = window.next() {
        if let Some(_args) = event.update_args() {
            grid.run_step(0.016); // Run one step (for ~60 FPS)
        }

        window.draw_2d(&event, |context, graphics, _device| {
            clear([0.1, 0.1, 0.1, 1.0], graphics); // Clear to dark gray

            // We'll create and call our drawing function here
            grid.draw_velocities(&context, graphics);
        });
    }

}
