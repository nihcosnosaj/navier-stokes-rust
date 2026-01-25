use std::mem;

pub struct FluidGrid {
    // Number of cells in x and y direction.
    pub nx: usize,
    pub ny: usize,

    // Cell width/height
    pub dx: f64,

    // Pressure (stored at cell centers)
    // Size: nx * ny
    p: Vec<f64>,
    p_prev: Vec<f64>, // Double buffer for pressure

    // Horizontal velocity (stored on vertical faces)
    // Size: (nx + 1) * ny
    pub u: Vec<f64>,
    u_prev: Vec<f64>, // Double buffer for u

    // Vertical velocity (stored on horizontal faces)
    // Size: nx * (ny + 1)
    pub v: Vec<f64>,
    v_prev: Vec<f64>, // Double buffer for v

    // Scratch space for divergence to avoid allocation in loop
    divergence: Vec<f64>,
}

impl FluidGrid {
    // A function to create a new, empty grid.
    pub fn new(nx: usize, ny: usize, dx: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            p: vec![0.0; nx * ny],
            p_prev: vec![0.0; nx * ny],
            u: vec![0.0; (nx + 1) * ny],
            u_prev: vec![0.0; (nx + 1) * ny],
            v: vec![0.0; nx * (ny + 1)],
            v_prev: vec![0.0; nx * (ny + 1)],
            divergence: vec![0.0; nx * ny],
        }
    }

    // A function to apply velocity decay.
    pub fn apply_velocity_decay(&mut self, decay: f64) {
        for u in &mut self.u {
            *u *= decay;
        }
        for v in &mut self.v {
            *v *= decay;
        }
    }

    // Convert 2D u-velocity index into 1D vector index.
    #[inline]
    pub fn u_idx(&self, i: usize, j: usize) -> usize {
        j * (self.nx + 1) + i
    }

    // Convert 2D v-velocity index to 1D vector index
    #[inline]
    pub fn v_idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    pub fn get_velocity(&self, x: f64, y: f64) -> (f64, f64) {
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
        // Copy current state to prev buffers to ensure boundaries are preserved
        // if we skip them in the loop.
        self.u_prev.copy_from_slice(&self.u);
        self.v_prev.copy_from_slice(&self.v);

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
                let idx = self.u_idx(i, j);
                self.u_prev[idx] = u_advected;
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
                let idx = self.v_idx(i, j);
                self.v_prev[idx] = v_advected;
            }
        }

        // Update the grid's velocity fields
        mem::swap(&mut self.u, &mut self.u_prev);
        mem::swap(&mut self.v, &mut self.v_prev);
    }

    #[inline]
    fn p_idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    // Step 2: Pressure Solve
    // Solves the Poisson equation to enforce incompressibility.
    fn solve_pressure(&mut self, _dt: f64) {
        // Implementation of the Jacobi iteration for the pressure solve.
        // This is where we calculate divergence and iterate to find p.
        let dx = self.dx;

        // Part 1: Calculate the divergence of the velocity field.
        // This is the right-hand side (RHS) of our Poisson equation.
        for j in 0..self.ny {
            for i in 0..self.nx {
                let u_right = self.u[self.u_idx(i + 1, j)];
                let u_left  = self.u[self.u_idx(i, j)];
                let v_top   = self.v[self.v_idx(i, j + 1)];
                let v_bot   = self.v[self.v_idx(i, j)];

                let d = (u_right - u_left + v_top - v_bot) / dx;
                let idx = self.p_idx(i, j);
                self.divergence[idx] = d;
            }
        }

        // Part 2: Iteratively solve for pressure using the Jacobi method.
        // We repeat this loop to let the pressure values settle.
        // Copy current p to p_prev to start with a good guess
        self.p_prev.copy_from_slice(&self.p);
        let num_iterations = 50; // More iterations = more accuracy
        for _ in 0..num_iterations {
            for j in 1..self.ny - 1 { // We only solve for interior pressure points
                for i in 1..self.nx - 1 {
                    let p_right = self.p[self.p_idx(i + 1, j)];
                    let p_left  = self.p[self.p_idx(i - 1, j)];
                    let p_top   = self.p[self.p_idx(i, j + 1)];
                    let p_bot   = self.p[self.p_idx(i, j - 1)];

                    let d = self.divergence[self.p_idx(i, j)];

                    // This is the discretized Poisson equation rearranged for p_i,j
                    let p_updated = (p_right + p_left + p_top + p_bot - d * dx * dx) / 4.0;
                    let idx = self.p_idx(i, j);
                    self.p_prev[idx] = p_updated;
                }
            }
            // Update the pressure field for the next iteration
            mem::swap(&mut self.p, &mut self.p_prev);
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

    fn diffuse_velocity(&mut self, viscosity: f64, dt: f64) {
        if viscosity == 0.0 {
            return;
        }
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        let a = viscosity * dt / (dx * dx);
        let iterations = 20;

        // Diffuse u
        // Initialize swap buffer with current state
        self.u_prev.copy_from_slice(&self.u);
        for _ in 0..iterations {
            for j in 1..ny - 1 {
                for i in 1..nx {
                    let idx = self.u_idx(i,j);
                    let idx_l = self.u_idx(i - 1, j);
                    let idx_r = self.u_idx(i+ 1, j);
                    let idx_b = self.u_idx(i, j - 1);
                    let idx_t = self.u_idx(i, j + 1);
                    // Jacobi iteration for diffusion
                    self.u_prev[idx] = (self.u[idx] + a * (self.u_prev[idx_l] + self.u_prev[idx_r] + self.u_prev[idx_b] + self.u_prev[idx_t]))
                        / (1.0 + 4.0 * a);
                }
            }
            // Swap buffers after each iteration so next iter reads updated values
            mem::swap(&mut self.u, &mut self.u_prev);
        }

        // Diffuse v
        self.v_prev.copy_from_slice(&self.v);
        for _ in 0..iterations {
            for j in 1..ny {
                for i in 1..nx - 1 {
                    let idx = self.v_idx(i, j);
                    let idx_l = self.v_idx(i - 1, j);
                    let idx_r = self.v_idx(i + 1, j);
                    let idx_b = self.v_idx(i, j - 1);
                    let idx_t = self.v_idx(i, j + 1);
                    self.v_prev[idx] = (self.v[idx] + a * (self.v_prev[idx_l] + self.v_prev[idx_r] + self.v_prev[idx_b] + self.v_prev[idx_t]))
                        / (1.0 + 4.0 * a);
                }
            }
            mem::swap(&mut self.v, &mut self.v_prev);
        }
    }

    pub fn run_step(&mut self, dt: f64, viscosity: f64) {
        self.diffuse_velocity(viscosity, dt);
        self.advect(dt);
        self.solve_pressure(dt);
        self.project(dt);
        self.set_boundaries();
        self.apply_velocity_decay(0.93);
    }

    pub fn add_velocity_at_pixel(&mut self, x: f64, y: f64, dx: f64, dy: f64) {
        // Convert window pixel coordinates to grid indices
        let i = (x / self.dx).floor() as isize;
        let j = (y / self.dx).floor() as isize;
        // Add velocity to a small area around the mouse
        for di in -1..=1 {
            for dj in -1..=1 {
                let ii = i + di;
                let jj = j + dj;
                if ii >= 0 && ii < self.nx as isize && jj >= 0 && jj < self.ny as isize {
                    let idx_v = self.v_idx(ii as usize, jj as usize);
                    let idx_u = self.u_idx(ii as usize, jj as usize);
                    // Scale the added velocity for effect
                    self.u[idx_u] += dx * 0.1;
                    self.v[idx_v] += dy * 0.1;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid_initialization() {
        let nx = 10;
        let ny = 10;
        let dx = 1.0;
        let grid = FluidGrid::new(nx, ny, dx);

        assert_eq!(grid.nx, nx);
        assert_eq!(grid.ny, ny);
        assert_eq!(grid.dx, dx);
        // Check buffer sizes
        // u is (nx + 1) * ny
        assert_eq!(grid.u.len(), (nx + 1) * ny);
        // v is nx * (ny + 1)
        assert_eq!(grid.v.len(), nx * (ny + 1));
    }

    #[test]
    fn test_indexing_logic() {
        let nx = 4;
        let ny = 4;
        let grid = FluidGrid::new(nx, ny, 1.0);

        // Test u_idx: index = j * (nx + 1) + i
        assert_eq!(grid.u_idx(0, 0), 0);
        assert_eq!(grid.u_idx(1, 0), 1);
        assert_eq!(grid.u_idx(0, 1), nx + 1);

        // Test v_idx: index = j * nx + i
        assert_eq!(grid.v_idx(0, 0), 0);
        assert_eq!(grid.v_idx(1, 0), 1);
        assert_eq!(grid.v_idx(0, 1), nx);
    }

    #[test]
    fn test_velocity_interpolation() {
        let mut grid = FluidGrid::new(10, 10, 1.0);
        // Set a uniform velocity field
        for u in &mut grid.u { *u = 10.0; }
        for v in &mut grid.v { *v = 20.0; }

        // Interpolation at any point should yield the uniform value
        let (u, v) = grid.get_velocity(5.5, 5.5);
        assert!((u - 10.0).abs() < 1e-6);
        assert!((v - 20.0).abs() < 1e-6);
    }

    #[test]
    fn test_add_velocity() {
        let mut grid = FluidGrid::new(10, 10, 1.0);
        
        // Add velocity at a specific pixel
        grid.add_velocity_at_pixel(5.0, 5.0, 100.0, 100.0);

        let (u, v) = grid.get_velocity(5.0, 5.0);
        assert!(u.abs() > 0.0 || v.abs() > 0.0, "Velocity should be non-zero after addition");
    }
}