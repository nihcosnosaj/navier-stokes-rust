use crate::grid::FluidGrid;
use piston_window::*;

pub struct FluidRenderer;

impl FluidRenderer {
    fn speed_to_color(speed: f64, max_speed: f64) -> [f32; 4] {
        // Map speed [0, max_speed] to blue->cyan->green->yellow->red gradient
        let t = (speed / max_speed).min(1.0).max(0.0);
        let (r,g,b) = if t < 0.25 {
            // Blue to Cyan
            (0.0, t * 4.0, 1.0)
        } else if t < 0.5 {
            // Cyan to Green
            (0.0, 1.0, 1.0 - (t - 0.25) * 4.0)
        } else if t < 0.75 {
            // Green to Yellow
            ((t - 0.5) * 4.0, 1.0, 0.0)
        } else {
            // Yellow to Red
            (1.0, 1.0 - (t - 0.75) * 4.0, 0.0)
        };
        [r as f32, g as f32, b as f32, 1.0]
    }

    /// Renders the density field as grayscale "smoke"
    pub fn draw_density(grid: &FluidGrid, c: &Context, g: &mut G2d) {
        // Find max density for normalization (prevents wash-out)
        let mut max_d = 0.01;
        for &d in &grid.density {
            if d > max_d { max_d = d; }
        }

        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let x = (i as f64) * grid.dx;
                let y = (j as f64) * grid.dx;
                
                // Fetch density from cell center
                let d = grid.density[j * grid.nx + i];
                
                // Map density to white/gray brightness
                // A common trick is to use (1.0 - exp(-d)) for a more "gaseous" look
                let val = (d / max_d).min(1.0) as f32;
                let color = [val, val, val, 1.0];

                rectangle(
                    color,
                    [x, y, grid.dx, grid.dx],
                    c.transform,
                    g,
                );
            }
        }
    }

    pub fn draw_velocities(grid: &FluidGrid, c: &Context, g: &mut G2d) {
        // Bright yellow, 3px thick
        // let line = line::Line::new([1.0, 1.0, 0.0, 1.0], 3.0);
        let scale = 2.0; // Larger scaling factor for arrow length

        // Find maximum speed for normalization.
        let mut max_speed = 1e-5;
        for j in 1..grid.ny - 1 {
            for i in 1..grid.nx - 1 {
                let x = (i as f64 + 0.5) * grid.dx;
                let y = (j as f64 + 0.5) * grid.dx;
                let (u, v) = grid.get_velocity(x, y);
                let speed = (u * u + v * v).sqrt();
                if speed > max_speed {
                    max_speed = speed;
                }   
            }
        }

        // Loop over a slightly smaller grid to avoid drawing on the boundaries
        for j in 1..grid.ny - 1 {
            for i in 1..grid.nx - 1 {
                let x = (i as f64 + 0.5) * grid.dx; 
                let y = (j as f64 + 0.5) * grid.dx;
                let (u, v) = grid.get_velocity(x, y);
                let speed = (u * u + v * v).sqrt();
                let color = Self::speed_to_color(speed, max_speed);
                let line = line::Line::new(color, 3.0);

                let x2 = x + u * scale;
                let y2 = y + v * scale;

                line.draw([x, y, x2, y2], &c.draw_state, c.transform, g);
            }
        }
    }

    pub fn draw_speed(grid: &FluidGrid, c: &Context, g: &mut G2d) {
        // 1. Find the maximum speed for normalization
        let mut max_speed = 1e-5;
        for j in 1..grid.ny - 1 {
            for i in 1..grid.nx - 1 {
                let x = (i as f64) * grid.dx;
                let y = (j as f64) * grid.dx;
                let (u, v) = grid.get_velocity(x + 0.5 * grid.dx, y + 0.5 * grid.dx);
                let speed = (u * u + v * v).sqrt();
                if speed > max_speed {
                    max_speed = speed;
                }
            }
        }

        // 2. Draw colored rectangles for each cell
        for j in 1..grid.ny - 1 {
            for i in 1..grid.nx - 1 {
                let x = (i as f64) * grid.dx;
                let y = (j as f64) * grid.dx;
                let (u, v) = grid.get_velocity(x + 0.5 * grid.dx, y + 0.5 * grid.dx);
                let speed = (u * u + v * v).sqrt();
                let color = Self::speed_to_color(speed, max_speed);

                rectangle(
                    color,
                    [x as f64, y as f64, grid.dx, grid.dx],
                    c.transform,
                    g,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_speed_to_color_gradient() {
        let max_speed = 100.0;

        // Test 0 speed (Blue)
        let c_zero = FluidRenderer::speed_to_color(0.0, max_speed);
        assert_eq!(c_zero, [0.0, 0.0, 1.0, 1.0]);

        // Test max speed (Red)
        let c_max = FluidRenderer::speed_to_color(max_speed, max_speed);
        assert_eq!(c_max, [1.0, 0.0, 0.0, 1.0]);

        // Test clamping (speed > max_speed should still be Red)
        let c_over = FluidRenderer::speed_to_color(max_speed * 2.0, max_speed);
        assert_eq!(c_over, [1.0, 0.0, 0.0, 1.0]);
    }
}